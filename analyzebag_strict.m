%% Used for all rosbag analysis. Used for experiment 3 of paper
% clear;
% close all; 

% Grasp values
CYLINDRICAL = 5; INDMID = 6;  PINCH = 7; ULNAR = 8; OPEN = 9; REST = 10; THUMBOPP = 11;
thisexpgrip = [THUMBOPP, CYLINDRICAL, INDMID, PINCH, ULNAR, OPEN];

% Objects are assigned numerical values
objectset = {'cylinderHandle::left_wheel', 'cross_joint_part::link', 'thin_rod'};
objectlabels = {'cylinder', 'cross joint', 'rectangle'};
objectnum = [0 1 2]; % taken from Nicolas' code
mapObj = containers.Map( objectset, objectnum);

NUMFINGERS = 4; % number of digits from the Allegrohand

% Tested conditions: either compliance or no compliance, grasp or release (no grasp) indicated to user
compliancegrasp_varnames = {'compnograsp', 'compgrasp', 'nocompgrasp', 'nocompnograsp'};

%Indices of the DoFs we care about - dictated by order of Nicolas' output
INDbase = 1; INDMCP = 2; INDPIP = 3; INDDIP = 4;
MIDbase = 5; MIDMCP = 6; MIDPIP = 7; MIDDIP = 8;
ENDbase = 9; ENDMCP = 10; ENDPIP = 11; ENDDIP = 12;
THUMBbase = 13; THUMBMCP = 14; THUMBPIP = 15; THUMBDIP = 16;

% Total contact achieved for each condition
catot = zeros(length(compliancegrasp_varnames), length(objectnum));
canumtot = zeros(length(compliancegrasp_varnames), length(objectnum));
rootdir = 'C:\Users\kzhuang\Documents\Visual Studio 2015\Projects\VMG30_Trentadue_collect\VMG30_Trentadue_collect\Data_out\';


%% Subject A1
% bags = {'201709071456','201709071538','201709071656','201709081129', '201709081702'  };
%% Subject A2
bags = {'201708021138','201708031635','201708031740','201708041003'};
%% Subject A3
% bags = {'201710311605'};
%% Subject B5
% bags = {'201712081503','201712081528','201801261147 '};
%% Subject B6
% bags = {'201712081855', '201801261410'}; 
%% Subject B7
% bags = {'201712141918'};  
 

%% Initialzing 
Cpercontact = [];
Cperconind = [];
Ctot = [];
Catot = [];
Htot = [];
grptot = [];
grp_betatot = [];
grp_perholdtot = [];
sesnum_tot = [];
objtot =  [];
hldtot =  [];
Envtot = [];

for bag = 1:length(bags)
    cd([rootdir bags{bag}]);
    files = dir;
    b={files.name};
    bagfiles = regexp(b,'bagdata\w*', 'match');
    for bagentry= 1:length(bagfiles)
        if ~isempty(bagfiles{bagentry})
            load(char(bagfiles{bagentry}));
        end
    end
    
    if ~exist('trial2', 'var')
        trial2 = trial;
    end
    
    %% Finding the different experimental conditions
    compliance = find(strcmp({trial2.compliance}, '3_per_finger'));
    nocompliance = find(strcmp({trial2.compliance}, 'no_contact'));
    grasp = find([trial2.instruction]);
    nograsp = find(~[trial2.instruction]);
    
    compnograsp = intersect(compliance, nograsp);
    compgrasp = intersect(compliance, grasp);
    nocompgrasp = intersect(nocompliance, grasp);
    nocompnograsp = intersect(nocompliance ,nograsp);
    
    %% Analyzing how long each grasp was held per trial as percentage of trial time
    camean = zeros(length(compliancegrasp_varnames), length(objectnum));
    
    % looping through each of the experimental conditions
    for c = 1:length(compliancegrasp_varnames)
        type = eval([compliancegrasp_varnames{c}]); % which experimental condition
        
        if ~isempty(type)
            contact = zeros(1, length(type));
            meanenv = zeros(1,length(type));
            chenv = zeros(size(trial2(1).EMGenvinterp, 1),length(type));
            conachieved = zeros(1, length(type));
            objectind = zeros(size(type));
            holdtot = [];
            holdgrp = [];
            
            contactottot = [];
            for i=1:length(type)            % looping through all trials of this condition
                object = trial(type(i)).object;             % which object grasped in words
                objectind(i) = mapObj(object);              % object index
                contacttot = trial(type(i)).totalcontacts;  % which pads were touching the object (value of 0,1, or 2 per pad)
                contacttot(contacttot>0)=1;                 % desired contacts = 2, contacts = 1, no contact = 0; 
                contactottot = [contactottot contacttot];   % all contacts for all trials of this type
                
                switch(objectind(i))
                    case 0 % cylindrical pad requirement for consideration as a successful grasp
                        tempgrasp = ( contacttot(INDPIP,:)>0 & contacttot(INDMCP,:)>0 &...
                            contacttot(MIDPIP,:)>0 & contacttot(MIDMCP,:)>0 &...
                            contacttot(ENDPIP,:)>0 & contacttot(ENDMCP,:)>0 &...
                            (contacttot(THUMBDIP,:)+ contacttot(THUMBPIP,:))>0);
                        grasp(i) = length(find(tempgrasp));
                    case 1 % cross-shaped pad requirement for consideration as a successful grasp
                        tempgrasp = (contacttot(INDDIP,:)>0 & contacttot(INDPIP,:)>0 &...
                            contacttot(MIDDIP,:)>0 & contacttot(MIDPIP,:)>0 &...
                            contacttot(ENDDIP,:)>0 & contacttot(ENDPIP,:)>0 &...
                            contacttot(THUMBDIP,:)>0)+ contacttot(THUMBPIP,:))>0);
                        grasp(i) = length(find(tempgrasp));
                    case 2 % thin bar pad requirement for consideration as a successful grasp
                        tempgrasp = (contacttot(INDDIP,:)>0 &contacttot(MIDDIP,:)>0 &...
                            (contacttot(THUMBDIP,:)+ contacttot(THUMBPIP,:))>0);
                        grasp(i) = length(find(tempgrasp));
                end
                
                contact(i) = grasp(i)/size(contacttot,2);           % fraction of trial time there was successful grasp     
                conachieved(i) = grasp(i)>0;                        % was there ever grasping (binary 1 or 0)
                
                % Calculating continuous hold time (uninterrupted)
                tempgrasp(find(diff(tempgrasp)==-1)+1)  = [];
                tempgrasp(find(diff(tempgrasp)==-1)+1)  = [];
                temp = num2str(tempgrasp,'%1d');
                holds = strsplit(temp, {'0', ' '});
                count = 1;
                holdtimes=[];
                for h = 1:length(holds)
                    if ~isempty(holds{h})
                        holdtimes(count) = length(holds{h});
                        count = count+1;
                    end
                end
                holdgrp = [holdgrp objectind(i)*ones(1, length(holdtimes))];    % indicates the object corresponding to the hold time
                holdtot = [holdtot holdtimes];                                  % hold time for the object of this trial
                
                emgenv = trial2(type(i)).EMGenvinterp;                          % EMG envelope for this trial
                meanenv(i) = mean(mean(emgenv(:,round(size(emgenv,2)/3):end))); % average for this trial is calculated from last 2/3 of trial
            end
            %     figure; boxplot(contact', 'Notch', 'on')
            eval(['ccon_' num2str(c) ' = contactottot;'])       % whether or not each pad made contact for all trials in this exp. condition
            eval(['c_' num2str(c) ' = contact;'])               % fraction of correct grasping for this experimental condition
            eval(['ca_' num2str(c) ' = conachieved;'])          % Whether or not contact was made for each trial in this exp. condition
            eval(['h_' num2str(c) '= holdtot;'])                % uninterrupted hold duration for each trial in this exp. condition
            eval(['env_' num2str(c) '= meanenv;'])              % Averaged emg enevelope during the trial

            % Data for each grasped object for this experimental condition
            for o = 1:length(objectnum)
                camean(c,o)= mean(conachieved(objectind ==(o-1)));                      % fraction of trials in which contact occurred
                catot(c,o) = catot(c,o)+sum(conachieved(objectind ==(o-1)));            % number of trials in which grasping occurred
                canumtot(c,o) = canumtot(c,o)+length(conachieved(objectind ==(o-1)));   % number of trials total
            end
            eval(['grp' num2str(c) ' = ' num2str(i) '*ones(size(contact));'])   % indicates number of trials for this experimental condition
            eval(['objgrp' num2str(c) ' = objectind;'])                         % indicates the object index for each trial
            eval(['hgrp' num2str(c) ' = holdgrp;'])                             % indicates object for each holdtime
        end
    end
    
    %% Data per experimental condition: 1 = release,COCO on; 2 = grasp, COCO on; 3 = grasp, COCO off, 4 = release, COCO off
    cpercon = [ccon_1 ccon_2 ccon_3 ccon_4]; 
    % Each trial is indexed by the experimental condition
    cpercon_ind = [ones(size(ccon_1,2),1); 2*ones(size(ccon_2,2),1);...
        3*ones(size(ccon_3,2),1); 4*ones(size(ccon_4,2),1)]';
    
    % Contact as percentage of trial time
    C=[c_1(:); c_2(:); c_3(:); c_4(:)];
    
    % Per trial binary contact data
    Ca=[ca_1(:); ca_2(:); ca_3(:); ca_4(:)];
    
    % Hold time per trial
    H =[h_1(:); h_2(:); h_3(:); h_4(:)];
    
    % Averaged EMG envelope per trial
    Env = [env_1(:); env_2(:); env_3(:); env_4(:)];
    
    % Keeps track of number of trials per experimental condition
    grp = [grp1(:); grp2(:); grp3(:); grp4(:)];
    
    % keeps track of object for each trial in each experimental condition
    obj = [objgrp1(:); objgrp2(:); objgrp3(:); objgrp4(:)];
    
    % Keeps track of object for each hold in each experimental condition
    hld = [hgrp1(:);hgrp2(:);hgrp3(:);hgrp4(:);];
    
    % Keeps track of fraction of time of correct grasp per trial per exp. condition
    grp_beta = [zeros(1,numel(c_1)),ones(1,numel(c_2)), 2*ones(1,numel(c_3)), 3*ones(1,numel(c_4))];
    
    % Keeps track of uninterrupted hold duration for each exp. condition
    grp_perhold = [zeros(1,numel(h_1)),ones(1,numel(h_2)), 2*ones(1,numel(h_3)), 3*ones(1,numel(h_4))];
    
    %% Aggregating data for all experimental conditions
    Cpercontact = [Cpercontact cpercon];
    Cperconind = [Cperconind cpercon_ind];
    Envtot= [Envtot; Env];
    Ctot = [Ctot; C];
    Catot = [Catot; C];
    Htot = [Htot; H];
    grptot = [grptot; grp];
    objtot =  [objtot; obj];
    hldtot = [hldtot; hld];
    sesnum_tot = [sesnum_tot bag*ones(size(grp_beta))];
    grp_betatot = [grp_betatot grp_beta];
    grp_perholdtot = [grp_perholdtot grp_perhold];
    
    numemgchannels = size(emgenv, 1);
    numtrialpts = 250;
    ComplianceEMG;
end

%% Single session plot
% figure; boxplot(C, {grp_beta, obj}, 'Notch', 'on')
% title('Percentage of Trial time there was categorical grasp');
% % figure; boxplot(C, {grp_beta, obj}, 'Notch', 'on')
% figure; bar(camean);
% figure; boxplot(H/30, {grp_perhold, hld}, 'Notch', 'on') % 30Hz sampling rate
% title('Hold times per grasp in seconds')
% % figure; boxplot(H, {grp_perhold, hld}, 'Notch', 'on')

%% Multisession plots

figure; boxplot(Ctot, {grp_betatot, objtot}, 'Notch', 'on')
title('Percentage of Trial time there was categorical grasp');

% Hypothesis testing for number of trial successes
for o = 1:3
    x = table([catot(2,o);catot(3,o) ],[canumtot(2,o)-catot(2,o);canumtot(3,o)-catot(3,o)],'VariableNames',{'Success','Failure'},'RowNames',{'COCO_ON','COCO_OFF'})
    [h,p,stats] = fishertest(x);
    [objectlabels(o) ' p-value with Fisher exact test for trial successes is: ' num2str(p)]
end

% Percentage of time of grasping success
figure; bar(catot./canumtot);
figure; plot( catot./canumtot, '.');
% Percentage of time of grasping success for only Grasping trials (no release)
figure; bar([(catot(3,:)./canumtot(3,:)); ((catot(2,:)./canumtot(2,:))-(catot(3,:)./canumtot(3,:)))]', 'stacked')
axis([0 4 0 1])

figure; boxplot(Htot/30, {grp_perholdtot, hldtot}, 'Notch', 'on') % 30Hz sampling rate
title('Hold times per grasp in seconds')

% figure; boxplot(H, {grp_perhold, hld}, 'Notch', 'on')
figure; boxplot(Envtot, {grp_betatot, objtot}, 'Notch', 'on')
title('Average EMG Amplitude during trial')

% looking at if there is any learning effect on number of trials in which grasping occurred
uniqueobj = unique(obj);
count =1;
uniquegrp_beta = [2,3];
figure;
hold on;

for uo= 1:length(uniqueobj)
    for ug = 1:length(uniquegrp_beta)
        uoinds = find(objtot == uniqueobj(uo));
        uginds = find(grp_betatot == uniquegrp_beta(ug));
        inds = intersect(uoinds, uginds);
        Ccell{count} = Ctot(inds);
        Ccellname{count} = ['Object ', objectset{uo}, ' Compliance ', num2str(uniquegrp_beta(ug))] ;
        plot(cumsum(Ctot(inds))./(1:length(inds))'); axis([-inf inf -2 4]);
        %             plot(cumsum(Ctot(inds))./repmat([1:length(inds)], )
        count = count+1;
    end
end
hold off;
legend(Ccellname)

%% Kruskal Wallis test on hold times
grpfind = find(grp_perholdtot == 1 | grp_perholdtot == 2);
for i = 2:length(unique(hldtot))
    objfind = find(hldtot ==(i-1));
    inds = intersect(grpfind, objfind);
    indscomp = intersect(find(grp_perholdtot ==1), objfind);
    indsnocomp = intersect(find(grp_perholdtot == 2), objfind);

    [P H STATS] = ranksum(Htot(indscomp)/30, Htot(indsnocomp)/30);
    P
end

figure; boxplot(Envtot, {grp_betatot}, 'Notch', 'on')
[P H STATS] = ranksum(Envtot(grp_betatot ==1), Envtot(grp_betatot ==2));
['EMG trials :P = ' num2str(P)]

%% Bar plot for number of contacts attained on average
tempc = find(Cperconind==2);
cond2 = Cpercontact(:,tempc);
tempc = find(Cperconind==3);
cond3 = Cpercontact(:,tempc);

tempc = [mean(cond3'); mean(cond2')];
figure; bar(tempc')
nonzer = find(mean(cond2')~=0);
figure; bar(tempc(:,nonzer)')
legend('Compliance OFF','Compliance ON')
axis([0 18 0 1])

%% Statistics for each contact pad (comparison between COCO on vs COCO off)
desiredcont = [2 3 4 6 7 8 10 11 12 15 16]; 
for o = 1:length(desiredcont)
    x = table([sum(cond2(desiredcont(o),:));sum(cond3(desiredcont(o),:))],[sum(cond2(desiredcont(o),:)==0);sum(cond3(desiredcont(o),:)==0)],'VariableNames',{'Success','Failure'},'RowNames',{'COCO_ON','COCO_OFF'});
%     y = table([sum(cond2(desiredcont(o),:));sum(cond3(desiredcont(o),:))],[length(cond2(desiredcont(o),:));length(cond3(desiredcont(o),:))],'VariableNames',{'Success','Failure'},'RowNames',{'COCO_ON','COCO_OFF'})
    [h,p,stats] = fishertest(x);
    [num2str(desiredcont(o)) ' p-value with Fisher exact test for trial successes is: ' num2str(p)]
end

