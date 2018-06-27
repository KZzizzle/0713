%% Does most of the data analysis for the Data collected in Experiment 2
% Obtains performance metrics and plots
% Must provide a cell array of experimental sessions for analysis. 


clear
close all

% Grasp types
THUMB =0; CYLINDRICAL = 5; INDMID = 6;  PINCH = 7; ULNAR = 8; OPEN = 9; REST = 10; THUMBOPP = 11;

randsucmattot =[];
ldasucmattot = [];
sucmattot = [];
sessiontot = [];
randsucmattot_tr = [];
sucmattot_tr = [];
trialsucmattot = [];
randtrialsucmattot = [];
realflextot = [];
predflextot= [];
rootdir = 'C:\Users\kzhuang\Documents\Visual Studio 2015\Projects\VMG30_Trentadue_collect\VMG30_Trentadue_collect\Data_out\';


%% Subject A1 using A2 params for fitting
% sessions = {'201709081215', '201709081232','201709081357', '201709081540'};
% sessions = {'201709081215', '201709081232','201709081540'};
%% Subject A1 
% sessions = {'201709061431', '201709061536', '201709070957','201709071032',...
%     '201709071204','201709071411','201709071430', '201709071604', '201709080954',...
%     '201709081103',  '201709081502' };
% sessions = {'201709081502' }; % for the confusion matrix
%% Subject A3 
% sessions = { '201710101037','201710101101', '201710311517', '201710311548'};
% sessions = { '201710311548'}; % for the confusion matrix
%% Subject A2 
% sessions = {'201708011310', '201708011529', '201708020956',...
%     '201708021010', '201708031538', '201708031602', '201708040919',...
%     '201708040943'};
% sessions = {'201708031538'};  % for the confusion matrix
%% Subject B6
% sessions = {'201712081734', '201801261330','201712081808','201801261352'}; 
% sessions = {'201712081808'};  % for the confusion matrix
%% Subject B5
% sessions = {'201712081416', '201712081439','201801261108','201801261124'}; 
% sessions = {'201801261124'}   % for the confusion matrix
%% Subject B3
 sessions = {'201712141751'} 
%% Subject B7
% sessions = {'201712141845','201712141906'}; 
% sessions = {'201712141906'}; % for the confusion matrix

r2  = zeros(length(sessions),12);
mse = zeros(length(sessions),12);

for s = 1:length(sessions)
    cd([rootdir sessions{s}]);
    files = dir;
    b={files.name};
    glovefiles = regexp(b,'vmg30\w*', 'match');
    firstflag = 0;
    
    kin_desc= {'Thumb MCP', 'Thumb PIP', 'Index MCP', 'Index PIP', 'Middle MCP',...
        'Middle PIP', 'Ring MCP', 'Ring PIP', 'Pinky MCP', 'Pinky PIP',...
        'Thumb Opp','Index Press', 'Middle Press', 'Ring Press', 'Pinky Press',...
        'Thumb-Index Abd/Add', 'Index-Middle Abd/Add', 'Middle-Ring Abd/Add', ...
        'Ring-Pinky Abd/Add', 'Thumb Cross', 'Palm Arch',...
        'Wrist Roll', 'Wrist Pitch', 'Wrist Yaw',...
        'Hand Roll', 'Hand Pitch', 'Hand Yaw',};
    
    %% Dividing into training and testing groups
    for i= 1:length(glovefiles)
        if (~isempty(glovefiles{i})&&firstflag==0)
            firstfile = (glovefiles{i});
            load(char(strcat(firstfile,'.txt')));
            realkin1 = eval(char(firstfile));
            trainedgrasps = unique([realkin1(:,end); REST]); % force it to count "rest" state as a grasp
            firstflag = 1;
        elseif (~isempty(glovefiles{i})&&firstflag~=0)
            secondfile = (glovefiles{i});
            load(char(strcat(secondfile,'.txt')));
            realkin2 = eval(char(secondfile));
            firstflag = 2;
        end
    end
    
    % if we have teleoperation-only sessions, set the vmg30 file to the testing set
    if firstflag ==1
        realkin2 =realkin1; 
    end
    realkin1(realkin1>1000)=1000;
    realkin2(realkin2>1000)=1000;
    realkin2(:, 1:end-1) = ((realkin2(:, 1:end-1))/...
        (max(max(realkin2(:,1:end-1))) - min(min(realkin2(:,1:end-1))))) ;
    
    %% Loading data files
    testingfile = regexp(b,'TESTING2\w*', 'match');
    for i= 1:length(testingfile)
        if (~isempty(testingfile{i}))
            load(char(strcat(testingfile{i},'.txt')));
            testing2 = eval(char(testingfile{i}));
        end
    end
    
    try
        LDAfile = regexp(b,'OUTPUTFeat_\w*', 'match');
        for i= 1:length(LDAfile)
            if (~isempty(LDAfile{i}))
                load(char(strcat(LDAfile{i},'.txt')));
                ldaf = eval(char(LDAfile{i}));
            end
        end
        ldapreds = ldaf(:,end-1:end);
        trainedinds = find(ismember(ldaf(:,end), trainedgrasps));
        ldapreds = ldaf(trainedinds,end-1:end);
        figure; plot(ldapreds);
        lda_suc(s) = length(find(ldapreds(:,1)==...
            ldapreds(:,2)))/size(ldapreds,1);
        
        lda_suc_trained(s) = length(find(ldapreds(trainedinds,1)==...
            ldapreds(trainedinds,2)))/length(trainedinds);
    catch
    end
    
    outputfile = regexp(b,'OUTPUTKinfilt\w*', 'match');
    for i= 1:length(outputfile)
        if (~isempty(outputfile{i}))
            load(char(strcat(outputfile{i},'.txt')));
            outputf = eval(char(outputfile{i}));
        end
    end
    
    % normalize the predictions
    outputf(:, 1:end-1) = ((outputf(:, 1:end-1))./...
        (repmat(((max(outputf(:,1:end-1))) - (min(outputf(:,1:end-1)))),...
        size(outputf, 1), 1))) ;
    
    %% Plot fit of the training data
    %     predtrain = testing(1:2:end,:);
    %     predtrainf = zeros(size(predtrain));
    %     for i = 1:size(predtrain,1); predtrainf(i,:) =smooth(predtrain(i,:), 20); end
    %     for i=1:size(predtrain,1)
    %         figure; plot(1:size(predtrain,2), [predtrainf(i,:); testing(i*2,:)]);
    %         axis([-inf inf -50 1050]);
    %         legend('Fitted Value', 'True Value');
    %         title(['Fitting Results With Training Data: ' kin_desc{i}]);
    %     end
    
    %% Calculate correlation coefficients of training data
    %     corrstr = zeros(1, size(predtrainf, 1));
    %     for i=1:size(predtrainf,1)
    %         correls = corrcoef(predtrainf(i,:), testing(i*2,:));
    %         corrstr(i) = correls(1,2);
    %     end
    
    %% Plot fit of data predicted in real time
    endind = min(size(realkin2,1),size(outputf,1));
    %         for i=1:size(realkin2,2)
    %             figure;
    %             plot(1:endind,[ outputf(1:endind,i) realkin2(1:endind,i)]');
    %             axis([-inf inf -50 1050]);
    %             title(['Real Time Running: ' kin_desc{i}]);
    %             legend('Fitted Value', 'True Value')
    %         end
    %
    %% Calculate correlation coefficients of real time data

    % occasionally, we tested a grasp that was not in the training set
    wasgrasptrained = find(ismember(realkin2(1:endind,end), trainedgrasps)); 
    for i=1:size(outputf,2)
        correls = corrcoef(realkin2(1:endind,i), outputf(1:endind,i));
        correls_trained = corrcoef(realkin2(wasgrasptrained,i), outputf(wasgrasptrained,i));
        corrs(s,i) = correls(1,2);
        corrs_tr(s,i) = correls_trained(1,2);
    end

    % Go through every other DoF since there are 2 per digit (MCP and PIP)
    for i=1:2:size(outputf,2)-1
        r_other = (realkin2(1:endind,i)+realkin2(1:endind,i+1));    % instructed kinematics
        o_other = (outputf(1:endind,i)+outputf(1:endind,i+1));      % predicted kinematics
        correls = corrcoef(r_other, o_other);
        corrs_other(s,i) = correls(1,2);
        mse_other(s,i) = mean(((o_other/(max(o_other)-min(o_other)))-(r_other/...
            (max(r_other)-min(r_other)))).^2);
        rms_other(s, i) = sqrt(mean(((o_other)-(r_other)).^2));
    end
    rsq = corrs.^2;
    
    %% Calculate correlation coefficients with SVM predictions
    svm_prelimscript; 

    %% Doing trial by trial analysis of decoding  performance
    NUMDOFS = 11;
    FS = 30;                    %30Hz
    grip = realkin2(:,end);             % last entry of kinematics data is the grasp type
    trials = find(diff(grip)~=0);       % is it the same trial
    trials = trials(diff(trials)>FS);   % make sure that the trial is at least a second long
    threshold = 0.17*(max(max(realkin2(:,1:end-1))) - min(min(realkin2(:,1:end-1)))); % 15 degrees within target angle
    
    trained = (ismember(grip, trainedgrasps));          % if the grasp was part of the trained grasps
    
    realflexed = zeros(length(trials), NUMDOFS);
    predflexed = zeros(length(trials), NUMDOFS);
    
    % Loop through all trials
    for i = length(trials)-1:-1:1
        
        trial(i).griptype = grip(trials(i)+2);          % actual grasp
        trial(i).wastrained = trained(trials(i)+2);     % was this grasp trained
        indices = trials(i) : trials(i+1)-1;            % number of indices in the trial
        trial(i).indices = indices;
        trial(i).length = length(indices);              % #time steps in this trial
        
        lda_threshold = .1;                             % must hold the grasp more than 10% of the trial
        intarget = zeros(1,NUMDOFS);                    % was the DoF within the instructed angle range
        intargettime = length(indices)*ones(1,NUMDOFS); % how long DoF was within the instructed angle range
        successes =[];
        for dof = 1:NUMDOFS
            realkin_hold = mean(realkin2(indices,dof)); % average DoF angle during trial
            realflexed(i, dof) = realkin_hold;          % instructed flexion amount for this DoF
            halftrial = round(length(indices)/4);       % analyze starting from this time index
            predflexed(i, dof) = mean(outputf(indices(halftrial:end),dof));  % predicted flexion amount for this DoF
            %----------------------------------------------------------------------
            % success if the difference between predicted and desired dof angle is below the threshold
            %----------------------------------------------------------------------
            success = find(abs(realkin_hold - outputf(indices,dof))<=threshold); % Find time points at which this DoF within specified angle limit
            if ~isempty(success)
                intarget(dof) = length(success)/length(indices);    % proportion of trial time inside target
                intargettime(dof) = success(1)/FS;                  % time at which success begins
                successes = [successes; success];
            else
                intarget(dof) =0;
            end
        end
        [j_success, ~] = hist(successes, unique(successes));            % number of joints within threshold range simultaneously
        [j_success, numjoints] = hist(j_success, unique(j_success));    % #time indices for each number of joints simultaneously within threshold
        sucinds = find(j_success<=(lda_threshold*length(indices)));     % time durations in which joint was not within requested bounds for long enough
        numjoints(sucinds) = [];            % if fewer than required amount of time held, remove those indices
        j_success(sucinds)=[];
        this_success = max(numjoints);      % maximum joints simultaneously correct for required amount of time

        if isempty(this_success); this_success = 0; end

        try
            trial(i).predgrip = g(ypred(indices)');    
            lda_success = length(find(g(ypred(indices)') == type(indices)))/length(indices);
            if lda_success <lda_threshold; lda_success =0; end
            trial(i).ldasuc = lda_success;
        catch
        end
        trial(i).trialsuc = this_success;           % maximum joints simultaneously correct for required amount of time
        trial(i).intarget = intarget;               % proportion of trial time inside target
        trial(i).intargettime = intargettime;       % time of beginning of correct positioning
    end
    
    %% LDA results (if they were calculated)
    try
        ldasucmat = [trial.lda_suc]';
    catch
        ldasucmat = zeros(size([trial.trialsuc]'));
    end
    
    sucmat      = reshape([trial.intarget], length(trial), length(trial(1).intarget));
    successes   = sum(reshape([trial.intarget], length(trial), length(trial(1).intarget)),2);
    
    %% Successful position metrics for only grasps that were in the training set
    trained_trials  = find([trial.wastrained]);
    sucmat_trained  = reshape([trial(trained_trials).intarget],...
        length(trained_trials), length(trial(1).intarget));                 
    successes_trained = sum(reshape([trial(trained_trials).intarget],...
        length(trained_trials), length(trial(1).intarget)),2);
    
    %% Randomized trials to assess success probability at random
    randtrials = randperm(length(trials)-1);    % randomize trial order
    for i = length(trials)-1:-1:1
        indices     = trials(i) : trials(i+1)-1;
        randindices = trials(randtrials(i)):trials(randtrials(i)+1);
        intarget = zeros(1,NUMDOFS);
        intargettime = length(indices)*ones(1,NUMDOFS);
        randsuccestot = [];
        for dof = 1:NUMDOFS
            if (dof<3)
                realkinposs = [ 0 0.35 1];
            elseif (dof<5)
                realkinposs = [0 0.7 1];
            else
                realkinposs = [0 1];
            end
            %             realkin_hold = mode(realkin2(randindices,dof));
            realkin_hold = realkinposs(randsample(length(realkinposs),1));
            randsuccess = find(abs(realkin_hold - outputf(indices,dof))<=threshold);
            %          if length(randsuccess)>(holdthreshold*length(indices))
            if ~isempty(randsuccess)
                intarget(dof) = length(randsuccess)/length(indices);
                intargettime(dof) = randsuccess(1)/FS;
            else
                intarget(dof) =0;
            end
            randsuccestot  =  [ randsuccestot; randsuccess];
        end
        [j_success, ~] = hist(randsuccestot, unique(randsuccestot));
        [j_success, numjoints] = hist(j_success, unique(j_success));
        randsucinds = find(j_success<=(lda_threshold*length(indices)));
        numjoints(randsucinds) = [];
        j_success(randsucinds) = [];
        this_randsuccess = max(numjoints);
        if isempty(this_randsuccess); this_randsuccess = 0; end;
        trial(i).randtrialsuc = this_randsuccess;
        trial(i).randintarget = intarget;
        trial(i).randintargettime = intargettime;
        
    end
    
    randsucmat = reshape([trial.randintarget], length(trial), length(trial(1).randintarget));
    randsuccesses = sum(reshape([trial.randintarget], length(trial), length(trial(1).randintarget)),2);
    
    randsucmat_trained = reshape([trial(trained_trials).randintarget], length(trained_trials), length(trial(1).randintarget));
    randsuccesses_trained = sum(reshape([trial(trained_trials).randintarget], length(trained_trials), length(trial(1).randintarget)),2);
    
    
    %% adding in results for each single session into larger matrix
    trialsucmattot = [trialsucmattot [trial.trialsuc]];
    randtrialsucmattot = [randtrialsucmattot [trial.randtrialsuc]];
    randsucmattot = [randsucmattot; randsucmat];
    sucmattot = [sucmattot; sucmat];
    sessiontot = [sessiontot s*ones(1, length(trials))];
    randsucmattot_tr = [randsucmattot_tr; randsucmat_trained];
    sucmattot_tr = [sucmattot_tr; sucmat_trained];
    
    realflextot = [realflextot; realflexed];
    predflextot = [predflextot; predflexed];
    
    
    %% Plotting fits
    figure;
    
    for k=1:NUMDOFS
        subplot(ceil(NUMDOFS/2), 2, k); plot((1:endind)/FS, [realkin2(1:endind,k), outputf(1:endind,k)]);
        axis([-inf/FS inf/FS -.5 2])
        r2(s,k)= 1-(sum((outputf(1:endind,k)-realkin2(1:endind,k)).^2)/sum((realkin2(1:endind,k)-mean(realkin2(1:endind,k))).^2));
        mse(s, k) = mean(((outputf(1:endind,k)/(max(outputf(1:endind,k))-min(outputf(1:endind,k))))-(realkin2(1:endind,k)/...
            (max(realkin2(1:endind,k))-min(realkin2(1:endind,k))) )).^2);
        rms(s, k) = sqrt(mean(((outputf(1:endind,k))-(realkin2(1:endind,k))).^2));
        %         mse(s,k) = sum((zscore(outputf(1:endind,k))-zscore(realkin2(1:endind,k))).^2)/length(1:endind);
    end
    
    %% State-ify the decoding into flexion/extension
    realstate = zeros(size(predflexed, 1), 11);
    predstate = zeros(size(predflexed, 1), 11);
    for i = 1:size(predflexed,2)
        if (mod(i,2)~=0)
            if (i~=size(predflexed,2))
                blah = realflexed(:,i) + realflexed(:,i+1);
            else
                blah = realflexed(:,i);
            end
            
            %% Flexion/Extension states for instructed movements
            real_flstate = zeros(size(realflexed(:,i)));
            real_exstate = zeros(size(realflexed(:,i)));
            real_flstate(blah >= (mean(blah)+(1*std(blah)))) = 1;    % makes threshold for flexion and categorizes DoF into 0 or 1
            real_exstate(blah <= (mean(blah)-(0.5*std(blah)))) = 1;  % makes threshold for extension and categorizes DoF into 0 or 1
            realstate(:,round(i/2)) = real_flstate;                  % one matrix for flexed or extended
            realstate(:,round(i/2)+round(size(realflexed,2)/2)) = real_exstate;                
            
            %% Flexion/Extension states for predicted movements
            pred_flstate = zeros(size(predflexed(:,i)));
            pred_exstate = zeros(size(predflexed(:,i)));
            if (i~=size(predflexed,2))
                blah = predflexed(:,i)+predflexed(:,i+1);
            else
                blah = predflexed(:,i);
            end
            pred_flstate(blah>= (1*mean(blah)+(0.5*std(blah)))) = 1;
            predstate(:,round(i/2)) = pred_flstate;                  % for flexion
            pred_exstate(blah<= (1.2*mean(blah)-(0*std(blah)))) = 1;
            predstate(:,round(i/2)+round(size(predflexed,2)/2)) = pred_exstate;                % for the extension
        end
    end
    
    %% Confusion matrix of flexion and extension
    diffs = predstate-realstate;
    for i = 1:size(diffs,2)
        inds = find(realstate(:,i)>0);
        confusion(s,:,i) = mean(diffs(inds,:));
    end 
end


%% Figures
% [p tbl stats] =anova1([sucmattot randsucmattot]);
% multcompare(stats);
% figure; boxplot([trialsucmattot(:)/NUMDOFS randtrialsucmattot(:)/NUMDOFS ldasucmattot], 'notch', 'on')
figure; boxplot([sucmattot(:) randsucmattot(:)], 'notch', 'on') %
figure; hist([sucmattot(:) randsucmattot(:)])

grp= ([zeros(size(sucmattot)) ones(size(randsucmattot))]);
grp1 = repmat([1:size(sucmattot,2)], size(sucmattot,1), 2);
grp2 = repmat([1:size(sucmattot,2)], size(sucmattot,1), 1);
figure;boxplot([sucmattot(:) randsucmattot(:)], {grp1(:), grp(:)}, 'notch', 'on')


%% Statistics on performance relative to random chance
%----------------------------------------------------------------------
% statistics for each dof (compared to random chance)
%----------------------------------------------------------------------
sucmattotvec = sucmattot(:);
randsucmatotvec = randsucmattot(:);
figure;
for i = 1:length(unique(grp2))
    doffind = find(grp2(:) ==(i));
    randsucmatotvec(doffind);
    [P,ANOVATAB,STATS] = kruskalwallis([sucmattotvec(doffind) randsucmattot(doffind)],[],  'off');
    P
    multcompare(STATS);
end
%----------------------------------------------------------------------
% statistics for all dofs together (compared to random chance and LDA)
%----------------------------------------------------------------------
% % [P,ANOVATAB,STATS] = kruskalwallis([randsample(trialsuc1dim, mintrials)/NUMDOFS,...
% %     randsample(randtrial1dim, mintrials)/NUMDOFS randsample(ldasucmattot, mintrials)],[],  'off');
% [P,ANOVATAB,STATS] = kruskalwallis([randsample(trialsuc1dim, mintrials)/NUMDOFS...
%     -median(randsample(randtrial1dim, mintrials)/NUMDOFS) ...
%     randsample(ldasucmattot, mintrials)-0.1],[],  'off');
% 
% P_alljoints = P
% figure; multcompare(STATS);


%% confusion matrix with both flexion and extension of each finger
confusionmlp1d = confusion(:);
confusionlda1d = confusion_lda(:); %must run svm_prelimscript

[P,ANOVATAB,STATS] = signrank(abs(confusionlda1d), abs(confusionmlp1d), 'alpha',0.01);
figure; boxplot((1-[abs(confusionlda1d) abs(confusionmlp1d)]).^2 , 'Symbol', '.')
title(P)
axis([0.75 2.25 -0.1 1.1])
% plot(cumsum(Ctot(inds))./(1:length(inds))')


