%% Caclulates average of EMG activity during grasping trials with Compliance control
% Called from analyzebag_strict.m - cannot be used standalone

%% For Reference....
% objectset = {'cylinderHandle::left_wheel', 'cross_joint_part::link', 'thin_rod'};
% compliancegrasp_varnames = {'compnograsp', 'compgrasp', 'nocompgrasp', 'nocompnograsp'};


objecttypes = {trial.object};   % object grasped during each trial
successful_grip = Ca;           % Whether or not grasp achieve for each trial

%% Loop through all experimental conditions and objects
figure;
for c = 1:length(compliancegrasp_varnames)
    compliances = eval(compliancegrasp_varnames{c});    % whether compliance control was on during the trials
    % loop through objects
    for o = 1:length(objectset)
        objects = find(ismember(objecttypes, objectset{o}));
        indices = (intersect(compliances, objects));    % indices of trials for desired exp condition and object
        numtrials = length(indices);                    % number of trials for desired exp cond. and object
        envelopes = zeros(numemgchannels, numtrialpts, numtrials);
        for t = 1:numtrials
            envelopes(:,:,t) = trial2(indices(t)).EMGenvinterp; % get EMG amplitudes from trial struct.
        end
        
        subplot(length(compliancegrasp_varnames), length(objectset), ((c-1)*(length(objectset)))+o);
%         imagesc((squeeze(mean(envelopes,3))')'); caxis([-1 1]); colormap  hot;    % heat map of EMG activity over trial time (one row per EMG channel)
        plot(mean(squeeze(mean(envelopes,3)))); axis([-inf inf -2 4]);              % line plot of traces of EMG activity over trial time
        title([num2str(length(indices)) ' trials averaged']);
    end
end

%%
% objecttypes = {trial.object};
% successful_grip = Ca;
% numemgchannels =7;
% numtrialpts = 250;
% % figure;
% for c = 1:length(compliancegrasp_varnames)
%
%     compliances = eval(compliancegrasp_varnames{c});
%
%     for o = 1:length(objectset)
%         objects = find(ismember(objecttypes, objectset{o}));
%         indices = intersect(compliances, objects);
%         envelopes =zeros(numemgchannels, numtrialpts, length(indices));
%
%         for t = 1:length(indices)
%             envelopes(:,:,t) = trial2(indices(t)).EMGenvinterp;
%             figure(t + ((bag-1)*10))
%             subplot(length(compliancegrasp_varnames), length(objectset), ((c-1)*(length(objectset)))+o);
%             imagesc(trial2(indices(t)).EMGenvinterp);caxis([-1 8]); colormap hot;
%         end
%
%         subplot(length(compliancegrasp_varnames), length(objectset), ((c-1)*(length(objectset)))+o);
%         imagesc((squeeze(mean(envelopes,3))')'); caxis([-1 8]); colormap hot;
% %         imagesc(trial2(indices(t)).EMGenvinterp);
%         title([num2str(length(indices)) ' trials averaged']);
%     end
% end


% figure; boxplot(Ctot, {grp_betatot, objtot}, 'Notch', 'on')