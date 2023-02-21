%% Targeting targeted memory reactivation: characteristics of cued reactivation in sleep
% Mahmoud E. A. Abdellahi, Anne C. M. Koopman, Matthias S. Treder & Penelope A. Lewis

% by Mahmoud Abdellahi
% emails: abdellahime@cardiff.ac.uk   m.eid@fci-cu.edu.eg
%%
function ids = lv_SO_Spindle_classification(features_test, ids, spindle)
% takes data formatted as in the feature extraction from headers and then a
% model is trained to predict correct and incorrect classification and then
% it is tested on unseen samples..

trained = 0;

if trained ==0
    %% prediction using slow oscillations
    path = 'D:\work\reactivation paper\sent to journal\data and code'; % data path ,example: 'D:\data and code'
    addpath([path '\data\generated from data and code'])
    if spindle ==0
        load Header_correct_exp Header_correct_exp; load Header_incorrect_exp  Header_incorrect_exp;
    else
        load Header_correct_exp2 Header_correct_exp2; load Header_incorrect_exp2  Header_incorrect_exp2;
        Header_correct_exp = Header_correct_exp2; Header_incorrect_exp = Header_incorrect_exp2;
    end

    sbj = [4 6 8 14 20 10 24 26	30 32 16 22]; % codes of participants
    labels=[];
    SO_features_correct = lv_SO_Spindle_FeatureExtraction(Header_correct_exp, spindle);
    SO_features_incorrect = lv_SO_Spindle_FeatureExtraction(Header_incorrect_exp, spindle);
    if spindle == 0
        feature_mat_correct=[]; feature_mat_incorrect=[];
        for i=1:numel(sbj)
            mx_dataClick_1 = SO_features_correct{i,1}; X = mx_dataClick_1{5, 10}.features;
            Z = X( ~isnan( mean(X,2) ) , : ); mx_dataClick_1{5, 10}.features = Z;
            feature_mat_correct{i,1} = [cos(mx_dataClick_1{5, 10}.features(:,1))  sin(mx_dataClick_1{5, 10}.features(:,1))   mx_dataClick_1{5, 10}.features(:,2:end)];

            mx_dataClick_1 = SO_features_incorrect{i,1}; X = mx_dataClick_1{5, 10}.features;
            Z = X( ~isnan( mean(X,2) ) , : ); mx_dataClick_1{5, 10}.features = Z;
            feature_mat_incorrect{i,1} = [cos(mx_dataClick_1{5, 10}.features(:,1))  sin(mx_dataClick_1{5, 10}.features(:,1))   mx_dataClick_1{5, 10}.features(:,2:end)];
        end
        features_correct = feature_mat_correct; features_incorrect = feature_mat_incorrect;
        for i=1:size(features_correct,1)
            labels{i,1} = [ones(size(cell2mat(features_correct(i,1)),1), 1) ; zeros(size(cell2mat(features_incorrect(i,1)),1), 1)];
            features{i,1} = [cell2mat(features_correct(i,1)) ; cell2mat(features_incorrect(i,1))];
        end
    else
        x = SO_features_correct'; x2 = SO_features_incorrect'; labels=[]; features=[];
        for i=1:size(x,1)
            feature_mat_correct{i,1} = cell2mat(x(i,:)); feature_mat_incorrect{i,1} = cell2mat(x2(i,:));
        end
        for i=1:size(feature_mat_correct,1)
            labels{i,1} = [ones(size(cell2mat(feature_mat_correct(i,1)),1), 1) ; zeros(size(cell2mat(feature_mat_incorrect(i,1)),1), 1)];
            features{i,1} = [cell2mat(feature_mat_correct(i,1)) ; cell2mat(feature_mat_incorrect(i,1))];
        end
    end
    %% classification
    % building train and test sets with all data
    featuresVals = cell2mat(features(: , 1)); clabels = cell2mat(labels(: , 1));
    TRAIN = [featuresVals clabels(:,1)];

    if spindle==0
        SO_threhsold = 75; vSOPeakTrough_pos = 5; % non SO rejection
        idx = TRAIN(:,vSOPeakTrough_pos)>SO_threhsold; TRAIN = TRAIN(idx,:);
    end
    GROUP_TRAIN = TRAIN(:,end);
    TRAIN(:,end) = []; % clearing the last column because it contains the labels
    [TRAIN,GROUP_TRAIN] = under_sample_data(TRAIN , GROUP_TRAIN); % balance classes

    % z-scoring
    TRAIN = zscore(TRAIN, [], 1);
    GROUP_TRAIN(GROUP_TRAIN==0) = 2;

    t = templateTree('NumVariablesToSample','all',...
        'PredictorSelection','interaction-curvature','Surrogate','all');
    rng(1);
    if spindle==0
        SO_prediction_mdl = fitcensemble(TRAIN,GROUP_TRAIN,'Method','Bag','NumLearningCycles',200,'Learners',t); % building the model
        save SO_prediction_mdl SO_prediction_mdl;
    else
        Spindle_prediction_mdl_new = fitcensemble(TRAIN,GROUP_TRAIN,'Method','Bag','NumLearningCycles',200,'Learners',t); % building the model
        save Spindle_prediction_mdl_new Spindle_prediction_mdl_new;
    end
else
    if spindle==0
        load SO_prediction_mdl SO_prediction_mdl; mdl = SO_prediction_mdl;
    else
        load Spindle_prediction_mdl Spindle_prediction_mdl; mdl = Spindle_prediction_mdl;
    end


    featuresVals2 = cell2mat(features_test(: , 1)); clabels2 = (ids(: , 1));
    TEST = [featuresVals2 clabels2(:,1)];
    if spindle==0
        SO_threhsold = 75; vSOPeakTrough_pos = 5; % non SO rejection
        idx = TEST(:,vSOPeakTrough_pos)>SO_threhsold; TEST = TEST(idx,:);
    end
    ids = TEST(:,end);
    TEST(:,end) = []; % clearing the last column because it contains the labels
    % z-scoring
    TEST = zscore(TEST, [], 1);

    [est,prob] = predict(mdl,TEST); % testing
    ids = ids(find(est==1));

end
%%
    function [TRAINhat , GROUP_TRAINhat] = under_sample_data(TRAIN,GROUP_TRAIN)
        TRAINc1 = TRAIN(GROUP_TRAIN==0 , :);
        TRAINc2 = TRAIN(GROUP_TRAIN==1 , :);
        if size(TRAINc1,1)>size(TRAINc2,1) % keeping near center
            k1 = size(TRAINc2,1);
            M = TRAINc1; C1 = mean(M ,1); D1= [];
            for k=1:size(M,1), D1(k,1) = (M(k,:)-C1)*(M(k,:)-C1)'; end
            [~,keepMrk1] = sort(D1,'ascend');
            TRAINc1 = TRAINc1(keepMrk1(1:k1),:);
        else
            k1 = size(TRAINc1,1);
            M = TRAINc2; C1 = mean(M ,1); D1= [];
            for k=1:size(M,1), D1(k,1) = (M(k,:)-C1)*(M(k,:)-C1)'; end
            [~,keepMrk1] = sort(D1,'ascend');
            TRAINc2 = TRAINc2(keepMrk1(1:k1),:);
        end
        TRAINhat = [TRAINc1 ; TRAINc2];
        GROUP_TRAINhat = [zeros(size(TRAINc1,1),1) ; ones(size(TRAINc2,1),1) ];
    end



end