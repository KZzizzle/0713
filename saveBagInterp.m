%% Mainly for synchronizing the bag file and the data output from Visual Studio
% This is necessary for plotting the EMG amplitude with respect to grasp
% and release trials in Gazebo


objectset = {'cylinderHandle::left_wheel', 'cross_joint_part::link', 'thin_rod'};
compliancegrasp_varnames = {'compnograsp', 'compgrasp', 'nocompgrasp', 'nocompnograsp'};

rootdir = 'C:\Users\kzhuang\Documents\Visual Studio 2015\Projects\VMG30_Trentadue_collect\VMG30_Trentadue_collect\Data_out\';
%% Subject A1
% bags = {'201709071456','201709071538','201709071656','201709081645', '201709081702'  };
%% Subject A2
% bags = {'201708021138','201708031635','201708031740','201708041003',};
%% Subject A3
% bags = {'201710311605'};
% bags = {'201708031740'};

for bag = 1:length(bags)
    cd([rootdir bags{bag}]);
    files = dir;
    b={files.name};
    bagfiles = regexp(b,'bagdata\w*', 'match');
    for bagentry= 1:length(bagfiles)
        if ~isempty(bagfiles{bagentry})
            loadfile = char(bagfiles{bagentry});
            load(loadfile);
        end
    end

    featurefile = regexp(b,'FEATURES\w*', 'match');
    for i= 1:length(featurefile)
        if (~isempty(featurefile{i}))
            load(char(strcat(featurefile{i},'.txt')));
            features = eval(char(featurefile{i}));
        end
    end
    
    
    outputfile = regexp(b,'OUTPUTKinfilt\w*', 'match');
    for i= 1:length(outputfile)
        if (~isempty(outputfile{i}))
            load(char(strcat(outputfile{i},'.txt')));
            prediction = eval(char(outputfile{i}));
        end
    end
    
    bagkindat = kinterp;
    bagkindat(isnan(bagkindat)) = 0;
    %% Number of EMG channels: may change per experiment
    numfeatures = 7;
    numemgchannels = size(features,1)/numfeatures;
    numtrialpts = 250;
    
    %% Hardcoded for the channels of VIE and virtual allegrohand
    bag2VIEmap = [16 2 3 6 7 10 11 13 15];
    VIE2bagmap = [1 3 4 5 6 7 8 11 2];
    bagdelay = 1;
    bagmapped_raw = bagkindat(bag2VIEmap,:)/pi*2000;
    predmapped = prediction(:,VIE2bagmap)';
    
    shorter_temp = min(length(sum(bagmapped_raw(:,bagdelay+1:end))),length(sum(predmapped)));
    shorter = min(shorter_temp, size(features, 2));
    bagmapped = bagmapped_raw(:,bagdelay:shorter+bagdelay-1);
    predmapped = predmapped(:,1:shorter);
    bagmappeds = zscore(sum(bagmapped(:,1:shorter)));
    predmappeds = zscore(sum(predmapped(:,1:shorter)));
    featurefiltered = zscore(envelope(features(1:numemgchannels, 1:shorter)', 10, 'peak'))';
    
%% Using findpeaks to get an idea of the lag between bag kinematics and VIE kinematics: must check by eye first...
    prom = .4;          % peak prominence
    smoothing = 150;    % smoothing window
    [predpks, ppki] = findpeaks(smooth(predmappeds, smoothing),'MinPeakProminence',prom, 'MinPeakDistance',200);
    [bagpks, bpki] =  findpeaks(smooth(bagmappeds,smoothing),'MinPeakProminence',prom,'MinPeakDistance',200);
    figure; findpeaks(smooth(predmappeds, smoothing),'MinPeakProminence',prom, 'MinPeakDistance',200);
    figure; findpeaks(smooth(bagmappeds,smoothing),'MinPeakProminence',prom,'MinPeakDistance',200);

    shorter2 = min(length(predpks), length(bagpks))-5;
    ind_diff = ppki(1:shorter2)-bpki(1:shorter2); 
    ind_diff_interp = interp1(ppki(1:shorter2), ind_diff,1:length(predmappeds),'linear', 'extrap');
    extraptime = (1:length(predmappeds)) - ind_diff_interp + bagdelay;
    prediction2 = zeros(size(prediction,2), length(globaltime));
    kinematics2 = bagmapped;
    figure;
    for j = 1:size(bagmapped,1)
        prediction2(j,:) = interp1(extraptime, predmapped(j,:), 1:length(globaltime));
        subplot(ceil(size(bagmapped,1)/2), 2, j); plot(1:shorter, [bagmapped_raw(j, 1:shorter); prediction2(j,1:shorter)]);
    end
    
    kinterp2 = bagmapped_raw;
    pinterp2 = pinterp; 
    
    emgenv = zeros(numemgchannels, length(globaltime));
    for j = 1:numemgchannels
        emgenv(j,:) = interp1(extraptime,featurefiltered(j,:),  1:length(globaltime));
    end
    emgenv(isnan(emgenv)) = 0;
    trial2 =trial;
    for i = 1:length(trial)
        [~, tindex_start] = min(abs(globaltime-trial(i).Tstart));
        [~, tindex_end] = min(abs(globaltime-trial(i).Tend));
        trial2(i).kinsbag = kinterp2(:, tindex_start:tindex_end);
        trial2(i).pressbag = pinterp2(:, tindex_start:tindex_end);
        startind = 31;
        numcontacts = 4;
        numfingers = 4;
        interpinds = linspace(tindex_start,tindex_end, numtrialpts);
        trial2(i).EMGenv = emgenv(:,tindex_start:tindex_end);
        emgenvinterp =zeros(numemgchannels, numtrialpts);
        for k = 1:numemgchannels
            emgenvinterp(k,:) = interp1( tindex_start:tindex_end,emgenv(k, tindex_start:tindex_end), interpinds);
        end
        emgenvinterp(isnan(emgenvinterp))=0;
        trial2(i).EMGenvinterp = emgenvinterp;
        trial2(i).totalcontactsbag = pinterp2(startind:end, tindex_start:tindex_end);
        trial2(i).totalpressuresbag = pinterp2(8:(8+16 -1), tindex_start:tindex_end);
        trial2(i).maxcontactbag = max(pinterp2(startind:end, tindex_start:tindex_end), 2);
        pnormal = pinterp2(startind:end, tindex_start:tindex_end);
        pnormal(pnormal>1) = 1;
        trial2(i).contacttimebag = sum(pnormal, 2);
        
        contactperfinger = zeros(4, length(tindex_start:tindex_end));
        for j=1:numfingers
            contactperfinger(j,:) = sum(pinterp2(startind+((j-1)*numcontacts):(startind+(j*numcontacts)-1), tindex_start:tindex_end));
        end
        trial2(i).contactperfingerbag = contactperfinger;
    end
    save(char([loadfile '_interp']), 'pinterp', 'cinterp', 'kinterp', 'linterp', 'iinterp', ...
        'pinterp2', 'kinterp2','prediction2',  'trial','trial2', 'globaltime', 'extraptime');
end