%% Parses the bag file from Nicolas' compliance controller in ROS simulation
% filename is the name of the .bag file you want to parse
% pinterp = tactile sensor pressures
% cinterp = compliance on/off (binary) 
% kinterp = kinematics decoded by MLP
% linterp = type of object grasped
% iinterp = instructed movement: grasp or release
% trial   = structure containing useful data for each grasp/release trial
% globaltime = time to which all other variables are interpolated

function [pinterp, cinterp, kinterp, linterp, iinterp, trial, globaltime] = parsebag(filename)

bagdata = rosbag(filename);

fs = 30;  % the update rate from Decoder

% possible desired contacts
gripset = {'full_tip', 'tripod', '3f', '2f', 'ulnar_pinch', '3_per_finger', 'no_contact', 'open'};
gripvalue = [0 1 2 3 4 11 12 12]; % taken from Nicolas' code
mapGrip = containers.Map(gripset,gripvalue);

% available objects for grasping
objectset = {'cylinderHandle::left_wheel', 'cross_joint_part::link', 'thin_rod'};
objectnum = [0 1 2]; % taken from Nicolas' code
mapObj = containers.Map(objectset,objectnum);


% Pressure sensors: contain both force and contact state data
press = select(bagdata, 'Topic', '/tactile_pressures');
pressdat = press.readMessages;
for i=1:size(pressdat) pdat(:,i) = pressdat{i}.Data; end;
ptime = press.timeseries.Time;
[C,IA,~] = unique(ptime);
ptime = C;
pdat = pdat(:,IA);


% Whether we are using the compliance controller or not. 3_per_finger = compliance on, no_contact = compliance off
compliance= select(bagdata, 'Topic', '/nl_command_parsed');
comptime = compliance.timeseries.Time;
compstring = compliance.readMessages;
cdat(1) = mapGrip('no_contact');
comptime = [bagdata.StartTime; comptime];
for i=1:size(compstring) cdat(i+1) = mapGrip(compstring{i}.Data);end;

% Kinematics of the robotic hand at the moment (from MLP sent over UDP)
kin = select(bagdata, 'Topic', '/desired_hand_state');
kindat = kin.readMessages;
ktime = kin.timeseries.Time;
for i=1:size(kindat) kdat(:,i) = kindat{i}.Position; end;
[C,IA,~] = unique(ktime);
ktime = C;
kdat = kdat(:,IA);

% Which link is being presented at the current time
link = select(bagdata, 'Topic', '/gazebo/set_link_state');
linkdat = link.readMessages;
linktime = link.timeseries.Time;
ldat = zeros(size(linktime));
for i=1:size(linkdat); ldat(i) = mapObj(linkdat{i}.LinkName); end;

% Whether to grasp or release the object at the current time: 1 = grasp, 0 = release
instruct = mod(1:length(ldat), 2);

% Experiment time
globaltime = bagdata.StartTime:1/fs:bagdata.EndTime;

% Interpolating all signals so we can have one standard time
pinterp = zeros(size(pdat, 1), length(globaltime));
for i=1:size(pdat,1); pinterp(i,:) = interp1(ptime, pdat(i,:), globaltime); end;
cinterp =interp1(comptime,cdat, globaltime, 'previous', 'extrap');
kinterp = zeros(size(kdat, 1), length(globaltime));
for i=1:size(kdat,1); kinterp(i,:) = interp1(ktime, kdat(i,:), globaltime); end;
linterp = interp1(linktime,ldat, globaltime, 'previous');
iinterp = interp1(linktime,instruct, globaltime, 'previous');

% parsing into trial structure
for i =  (length(ldat)-1):-1:1
    trial(i).object = objectset{ldat(i)+1};
    
    trial(i).orientation    = linkdat{i}.Pose.Orientation;  % orientation of the object
    trial(i).Tstart         = linktime(i);                  % start time of trial
    trial(i).Tend           = linktime(i+1);                % end time of trial
    [~, tindex_start]       = min(abs(globaltime-linktime(i)));     % index of start in globaltime
    [~, tindex_end]         = min(abs(globaltime-linktime(i+1)));   % index of end 
    trial(i).time           = globaltime(tindex_start:tindex_end);  % time points of trial
    trial(i).kins           = kinterp(:, tindex_start:tindex_end);  % robot hand kinematics during trial
    trial(i).press          = pinterp(:, tindex_start:tindex_end);  % pressure sensation during trial
    
    startind    = 31;
    numcontacts = 4;
    numfingers  = 4; 
    trial(i).totalcontacts  = pinterp(startind:end, tindex_start:tindex_end);   % binary whether or not there is contact on a pad 
    trial(i).totalpressures = pinterp(8:(8+16 -1), tindex_start:tindex_end);    % 16 dof for hand, 8 for arm. Nicolas' code starts at 0 for indexing instead of 1 for Matlab
    trial(i).maxcontact     = max(pinterp(startind:end, tindex_start:tindex_end), 2);   % maximum number of contacts within a trial
    pnormal                 = pinterp(startind:end, tindex_start:tindex_end);
    pnormal(pnormal>1)      = 1;
    trial(i).contacttime    = sum(pnormal, 2);                                  % duration of contact 
    
    % how many contacts on a single fingure during the trial
    contactperfinger = zeros(4, length(tindex_start:tindex_end));
    for j=1:numfingers
         contactperfinger(j,:) = sum(pinterp(startind+((j-1)*numcontacts):(startind+(j*numcontacts)-1), tindex_start:tindex_end));
    end
    trial(i).contactperfinger = contactperfinger;
    
    % list of desired contacts (skip the most medial one for all trials) - one set of 4 values for each finger
    if (iinterp(tindex_start+1)==1 && ldat(i)~=2)
        trial(i).desired_contacts = ...
            repmat([0 1 1 1 0 1 1 1 0 1 1 1 0 0 1 1]', 1, length(tindex_start:tindex_end));
    elseif (iinterp(tindex_start+1)==1 && ldat(i)==2)
        trial(i).desired_contacts = ...
            repmat([0 1 1 1 0 1 1 1 0 0 0 0 0 0 1 1]', 1, length(tindex_start:tindex_end));
    else
        trial(i).desired_contacts = ...
            zeros(16, length(tindex_start:tindex_end));
    end

    % How much time in a trial before contact starts
    contactperfinger(contactperfinger>1)=1;
    delay_ind = find(sum(contactperfinger)>2);
    if ~isempty(delay_ind)
        trial(i).contactdelay = delay_ind(1)/fs;
    else
        trial(i).contactdelay = -1;
    end
    
    trial(i).compliance = gripset{gripvalue == cinterp(tindex_start+1)};    % compliance on or off
    trial(i).instruction = iinterp(tindex_start+1);                         % instruction was grasp or release
end
%
% map_emg_to_allegro_3(emg_joints):
% allegro_joints = np.zeros(16)
% # index
% allegro_joints[0] = 0.0 % adduction/abduction
% allegro_joints[1] = emg_joints[1]
% allegro_joints[2] = emg_joints[2]
% allegro_joints[3] = emg_joints[3]
% # middle: get it from the index
% allegro_joints[4] = 0.0 % adduction/abduction
% allegro_joints[5] = emg_joints[5]
% allegro_joints[6] = emg_joints[6]
% allegro_joints[7] = emg_joints[7]
% # pinky
% allegro_joints[8] = 0.0  % adduction/abduction
% allegro_joints[9] = emg_joints[13]
% allegro_joints[10] = emg_joints[14]
% allegro_joints[11] = emg_joints[15]
% # thumb
% allegro_joints[12] = emg_joints[16]
% allegro_joints[13] = -0.1
% allegro_joints[14] = emg_joints[18]
% allegro_joints[15] = emg_joints[19]
