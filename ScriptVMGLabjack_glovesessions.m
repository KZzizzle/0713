%% Does "confusion matrix" analysis and accuracy of decoding with the kinematic glove.
% Used for Figure 3 

close all

THUMB =0; CYLINDRICAL = 5; INDMID = 6;  PINCH = 7; ULNAR = 8; OPEN = 9; REST = 10; THUMBOPP = 11;
rootdir = 'C:\Users\kzhuang\Documents\Visual Studio 2015\Projects\VMG30_Trentadue_collect\VMG30_Trentadue_collect\Data_out\';

sessions = {'201702091759'}; % Subject B3


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
    
    %% Loading data
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
    if firstflag ==1
        realkin2 =realkin1; % if we have teleoperation-only sessions, set the vmg30 file to the testing set
    end
    realkin1(realkin1>1000)=1000;
    realkin2(realkin2>1000)=1000;
    
    %% Loading more data
    testingfile = regexp(b,'TESTING2\w*', 'match');
    for i= 1:length(testingfile);
        if (~isempty(testingfile{i}))
            load(char(strcat(testingfile{i},'.txt')));
            testing2 = eval(char(testingfile{i}));
        end
    end

    outputfile = regexp(b,'OUTPUTKinfilt\w*', 'match');
    for i= 1:length(outputfile)
        if (~isempty(outputfile{i}))
            load(char(strcat(outputfile{i},'.txt')));
            outputf = eval(char(outputfile{i}));
        end
    end
    
    outputf     = outputf(1:end-round(size(outputf, 1)*.15), :);
    realkin2    = realkin2(1:size(outputf, 1), :);
    
    realstate = zeros(size(outputf, 1), 6);
    predstate = zeros(size(outputf, 1), 6);

    figure(51);
    
    %% Looping through all DoFs
    for i=1:12
        subplot(6, 2, i); plot((1:size(realkin2,1))/30, [realkin2(:,i+1) outputf(:,i)]);
        axis([-inf inf -50 2050]);
        legend( 'True Value','Fitted Value');
        title(['Fitting Results With Real Time Data: ' kin_desc{i}]);
        
        %================ state-ify the decoding into flexion/extension ==========================
        if (mod(i,2)~=0)
            real_flstate = zeros(size(realkin2(:,i+1)));
            real_exstate = zeros(size(realkin2(:,i+1)));
            blah = realkin2(:,i+1)+realkin2(:,i+2);
            real_flstate(blah>= (mean(blah)+(1.7*std(blah)))) = 1;
            real_exstate(blah<= (mean(blah)-(0*std(blah)))) = 1;
            realstate(:,round(i/2)) = real_flstate;                     % for the flexion
            realstate(:,round(i/2)+6) = real_exstate;                   % for the extension
            
            pred_flstate = zeros(size(outputf(:,i)));
            pred_exstate = zeros(size(outputf(:,i)));
            blah = outputf(:,i)+outputf(:,i+1);
            pred_flstate(blah>= (mean(blah)+(.5*std(blah)))) = 1;       % predicted flexion
            predstate(:,round(i/2)) = pred_flstate;   
            pred_exstate(blah<= (1.1*mean(blah)-(0*std(blah)))) = 1;    % predicted extension
            predstate(:,round(i/2)+6) = pred_exstate;                  
        end
        %==================================================================
        
        %% Performance Metrics
        correls = corrcoef(realkin2(:,i+1)+realkin2(:,i+2), outputf(:,i)+outputf(:,i+1));
        corrs_other(s,i) = correls(1,2);
        %         r2(s,i)= 1-(sum((outputf(:,i)+outputf(:,i+1)-realkin2(:,i+1)-realkin2(:,i+2)).^2)/...
        %             sum((realkin2(:,i+1)+realkin2(:,i+2)-mean(realkin2(:,i+1)+realkin2(:,i+2))).^2));
        r2(s,i)  = 1-(sum((outputf(:,i)-realkin2(:,i+1)).^2)/...
            sum((realkin2(:,i+1)-mean(realkin2(:,i+1))).^2));
        correls_other = corrcoef(realkin2(:,i+1), outputf(:,i));    % correlation with every other DoF
        corrs(s,i) = correls_other(1,2);
        nmse(s, i) = mean(((outputf(:,i)/(max(outputf(:,i))-min(outputf(:,i))))-(realkin2(:,i+1)/...
            (max(realkin2(:,i+1))-min(realkin2(:,i+1))) )).^2);
%         nmse(s, i) = mean(((zscore(outputf(:,i))- zscore(realkin2(:,i+1))).^2)/(max(zscore(outputf(:,i)))-min(zscore(outputf(:,i)))));
%         r2(s,i)= 1-(sum((outputf(:,i)-realkin2(:,i+1)).^2)/sum((realkin2(:,i+1)-mean(realkin2(:,i+1))).^2));
    end
    %% Confusion matrix for flexion/extension of each DoF
    diffs = predstate-realstate;
    for i = 1:size(diffs,2)
        inds = find(realstate(:,i)>0);
        confusion(s,:,i) = mean(diffs(inds,:));
    end
    
end
confusion(isnan(confusion))=0;
figure; imagesc(squeeze(mean(confusion, 1)))