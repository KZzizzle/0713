%% Prediction of joint angles with SVM or LDA
% Implements an algorithm using the features that the MLP uses and also is
% trained on the same kinematic data as MLP

%% Format data into the right dimensions for training and testing
featurefile = regexp(b,'FEATURES\w*', 'match');
for i= 1:length(featurefile)
    if (~isempty(featurefile{i}))
        load(char(strcat(featurefile{i},'.txt')));
        features = eval(char(featurefile{i}));
    end
end

%% Some files are zero-padded in the buffer so try to find the end of the session
endpoint = strfind(sum(realkin1, 2)', zeros(1, round(0.1*size(realkin1,1))));
if isempty(endpoint); endind1 = size(realkin1, 1);
else; endind1 = endpoint(1); end

endind2 = min(size(realkin2,1),size(outputf,1));
feature1 = features(1:round(size(features,1)/2), 1:endind1);
feature2 = features(round(size(features,1)/2)+1:end, 1:endind2);

%% Load data
kintrain1 = realkin1(1:endind1,:)';
kintrain2 = realkin2(1:endind2,:)';

trlen= size(feature1,2);
tslen= size(feature2,2);

if firstflag~=1
    FeaturesTr = feature1';
    FeaturesTs = feature2';
else
    trlen= round(size(feature2,2)*0.6);
    tslen= size(feature2,2)-trlen-1;
    
    FeaturesTr = feature2(:, 1:trlen)';
    FeaturesTs = feature2(:, end-tslen+1:end)';
    kintrain1 = realkin1(1:trlen,:)';
    kintrain2 = realkin2(end-tslen+1:end,:)';
end

%% The actual fit and prediction part
% svmpred = zeros(size(kintrain2));
% for i = 1:size(kintrain1, 1)
%     MDL = fitlm( FeaturesTr, kintrain1(i,:));
%     YP = predict(MDL, FeaturesTs);
%     YP = smooth(YP, 10);
%     svmpred(i, :) = YP;
% %     figure; plot([YP 1000*kintrain2(i,:)'])
% end

%% LDA prediction per joint
% stateout = zeros(endind, 6);
kintr = zeros( 6, trlen);
kints = zeros( 6 ,tslen);
count =1;
%% Loop through the DoFs
for i=1:2:size(kintrain1,1)-2
    %         stateout(:,count) = (outputf(1:endind,i)+outputf(1:endind,i+1))/2;
    kintr(count,:) = (kintrain1(i,:)+kintrain1(i+1,:))/2000;
    kints(count,:) = (kintrain2(i,:)+kintrain2(i+1,:))/2;
    count = count+1;
end
kintr(count,:) = kintrain1(end-1,:)/1000;
kints(count,:) = kintrain2(end-1, :)/1;
clear ('statetr_fl', 'statetr_ex', 'statets_fl', 'statets_ex')
%% Stateify the kinematics into flexion and extension
for i = 1:size(kintr,1)
    statetr_fl(i,:) = kintr(i,:) >= (mean(kintr(i,:))+0.5*std(kintr(i,:)));
    statetr_ex(i,:) = kintr(i,:) <= (1.2*mean(kintr(i,:))-(0*std(kintr(i,:))));
    statets_fl(i,:) = kints(i,:) >=(mean(kints(i,:))+0.5*std(kints(i,:)));
    statets_ex(i,:) = kints(i,:) <= (1.2*mean(kints(i,:))-(0*std(kints(i,:))));
end
% statetr_fl = kintr>(mean(kintr)+0.5*std(kintr));
% statetr_ex = kintr<0.4;
% statets_fl = round(kints);
% statets_ex = kints<0.4;

state_tr = [statetr_fl; statetr_ex];
state_ts = [statets_fl; statets_ex];

%% Fit LDA
state_pred = zeros(size(kintrain2));
for i = 1:size(kintrain1, 1)
    MDL = fitcdiscr( FeaturesTr, state_tr(i,:));
    YP = predict(MDL, FeaturesTs);
    state_pred(i, :) = YP;
end

%% LDA confusion matrix
diffs = state_ts-state_pred;
for i = 1:size(diffs,1)
    inds = find(state_ts(i,:)>0);
    confusion_lda(s,:,i) = mean(diffs(:,inds)');
end

grptr =([zeros(size(sucmattot_tr)) ones(size(randsucmattot_tr))]);
grp1tr = repmat(1:size(sucmattot_tr,2), size(sucmattot_tr,1), 2);
mean(mean(abs(confusion_lda(s,:,:))))
caxis([-1 1])
