function [ txtAuc ] = lv_classify(cfg)
% takes trn and tst data in case on timextime and returns the
% txtPerformance .. if one set in provided then it's either timextime or
% axtime according to cfg.method and it will be CV on the same dataset ..
%EXAMPLE: txt:
% cfg=[];
% cfg.method = 'timextime';
% cfg.classifier_type = 'lda';
% cfg.trn = sleep; %has labels in trialinfo
% cfg.tst = im; cfg.folds=nan;
% [ txtAuc ] = lv_classify(cfg);

if isfield(cfg,'weights'), observation_weights = cfg.weights; else observation_weights = ones(size(cfg.trn.trial,1),1); end
method = {'axtime','timextime'};
classifier_type = {'lda','svm','random_forest','naive_bayes'};
if ~isfield(cfg,'method')
    idx = listdlg('ListString',method, 'ListSize',[450,450]);
    idx2 = listdlg('ListString',classifier_type, 'ListSize',[450,450]);
    cfg.method = char( method(idx) );
    cfg.classifier_type = char( classifier_type(idx2) );
end
fprintf(['\n Performing: ' cfg.method ' classification \n']);

if ~isfield(cfg,'perf_measure'), cfg.perf_measure = 'auc'; end

TRAIN = cfg.trn.trial;
GROUP_TRAIN = cfg.trn.trialinfo(:); % always column vector containing labels

if isfield(cfg,'every_channel'), channels = num2cell(1:size(TRAIN,2)); else % if classification is on each channel then one feature will be used
    channels = mat2cell( 1:size(TRAIN,2),1, size(TRAIN,2) ); end

if isnan(cfg.folds) % if folds is set with a number then we are dealing with the same set trn and we have cross validation
    GROUP_TEST = cfg.tst.trialinfo(:);
    TEST  = cfg.tst.trial; folds=nan; else
    folds = cfg.folds;
end

if cfg.do_parallel==1, p = gcp; workers=p.NumWorkers; else, workers=0; end

% takes 3D of training and testing sets and returns timextime AUC
% tip: elle 3ayz tefredo zy el time 7oto fl dim el tany w talt w elle 3ayzo yfdal zy el features 5aleh fl awal ...

if folds==0 % folds=0: LeaveOut, nan:two datasets, number: no. folds on the same dataset
    rng(1), CVO = cvpartition(size(TRAIN,1),'LeaveOut');
    warning('Doing z-scoring with stats of big fold because leave one out, so we only have one test trial..');
elseif isnan(folds)
    CVO.NumTestSets=1;
else
    rng(1)      % CV on the same dataset
    CVO = cvpartition(GROUP_TRAIN,'k',folds);
end


TRAIN_hold=TRAIN; if isnan(cfg.folds), TEST_hold=TEST; end
for fi = 1:CVO.NumTestSets % folds
    fprintf(['\n Fold: ' num2str(fi) '\n']);
    for ch=1:length(channels) % possible seperate channel classification
        TRAIN=TRAIN_hold(:,cell2mat(channels(ch)),:); if isnan(cfg.folds), TEST=TEST_hold(:,cell2mat(channels(ch)),:); end
        
        if isnan(folds) % txt for two datasets
            fold_TRAIN = TRAIN; fold_TEST = TEST;
            fold_TRAIN_GROUP = GROUP_TRAIN;
            fold_TEST_GROUP = GROUP_TEST;
        else
            rng(1)      % CV on the same dataset
            fold_TRAIN = TRAIN(CVO.training(fi),:,:);
            fold_TEST = TRAIN(CVO.test(fi),:,:); % really stratified -> sum(GROUP(CVO.test(i))==1)
            fold_TRAIN_GROUP = GROUP_TRAIN(CVO.training(fi));
            fold_TEST_GROUP = GROUP_TRAIN(CVO.test(fi));
            % z-scoring the fold
            if size(fold_TEST,1)>1,fold_TRAIN = zscore(fold_TRAIN, [], 1); % we zscore if folds not leave one out
                fold_TEST = zscore(fold_TEST, [], 1); else % leave one out is zscored with mu and sigma of the big fold
                fold_TEST = (fold_TEST-mean(fold_TRAIN,1))./std(fold_TRAIN,[],1);
                fold_TRAIN = zscore(fold_TRAIN, [], 1);
            end
            lv_progress(ch, length(channels), 'channels:');
        end
        if strcmp(cfg.method,'timextime')
            % forming a trltime dimension for faster classification
            dim = size(fold_TEST);
            temp = permute(fold_TEST,[2 1 3]);
            temp = reshape(temp, dim(2),[]); % ch x trltime
            fold_TEST = temp'; %  trltime x ch
            dimension = ones(1,size(fold_TRAIN,3)); % just not to make a condition inside the parfor to be faster
            % this will make code always access the first third dim. and that's
            % what we want because the data is trltime x ch
            if length(dim)==2, dim(3)=1; end % if data is 2d assume we have one time pt in the third dim to handle 2d
        elseif strcmp(cfg.method,'deep') % deep is special because it takes the whole 3d data and classify it
            % as a sequence in time
            fold_TRAIN_GROUP_temp = fold_TRAIN_GROUP; fold_TRAIN_GROUP=[];
            if strcmp(cfg.classifier_type(1),'TCN')==0
                classes_labels = [repmat( categorical(cellstr('c1')),1,size(fold_TRAIN,3));repmat( categorical(cellstr('c2')),1,size(fold_TRAIN,3));...
                    repmat( categorical(cellstr('c3')),1,size(fold_TRAIN,3));repmat( categorical(cellstr('c4')),1,size(fold_TRAIN,3))];
                for i=1:size(fold_TRAIN,1), temp{1,i} = squeeze(fold_TRAIN(i,:,:));
                    fold_TRAIN_GROUP{1,i} = classes_labels(fold_TRAIN_GROUP_temp(i),:); end
                if strcmp(cfg.classifier_type(1),'1d_cnn')==1, fold_TRAIN_GROUP=categorical(fold_TRAIN_GROUP_temp); end
            else 
                for i=1:size(fold_TRAIN,1), temp{1,i} = squeeze(fold_TRAIN(i,:,:));
                    fold_TRAIN_GROUP{1,i} = categorical(repmat((fold_TRAIN_GROUP_temp(i)),1,size(fold_TRAIN(i,:,:),3))); end
            end
            fold_TRAIN = temp; clear temp;


            dim=size(fold_TEST);
            for i=1:size(fold_TEST,1), temp{1,i} = squeeze(fold_TEST(i,:,:)); end
            fold_TEST = temp; clear temp;
            cfg.perf_measure = 'acc';
            dimension = 1;
        else
            dimension = 1:size(fold_TRAIN,3); dim = size(fold_TEST); dim(3)=1; % to make it one time pt because we don't test many time pts just the corresponding from test
            % such that we access the same time point from trn and tst and
            % that's what we need because the data is here still 3d

        end
        
        if isfield(cfg,'weights'), observation_weights = cfg.weights; else observation_weights = ones(size(fold_TRAIN,1),1); end % weights are relative to the training fold
        % datatypes conversion
        if strcmp(class(fold_TRAIN),'tall'), fold_TRAIN=gather(fold_TRAIN); TRAIN=gather(TRAIN); fold_TEST=gather(fold_TEST); end
        if numel(fold_TRAIN)>(2154*13*1601), fold_TRAIN=single(fold_TRAIN); TRAIN=single(TRAIN); fold_TEST=single(fold_TEST); end % size limit (2154*13*1601) if exceeded will convert to single

        if length(cfg.classifier_type)==3, if strcmp(cfg.classifier_type{2},'domain_align'), cfg.classifier_type{3}=dim; end, end
        % because in domain align all timepts of train are aggregated but we need to separate them later on to adapt to trn's domain
        if strcmp(cfg.perf_measure,'dval')==1, fold_dval = CVO.test(fi); else fold_dval = []; end % this one was if fi>1, fold_dval = CVO.test(fi); else fold_dval=1; end but it was giving error when we use 10fold to get dvals
        %         progressbar = ParforProgressbar(size(fold_TRAIN,3),'title', 'Classification progress');
%         parfor (j=1:size(fold_TRAIN,3), workers) % no. workers,, will be 0 if not parallel and max available if parallel
                                            for j=1:size(fold_TRAIN,3)

            [outclass,~,hh] = low_lvl_classify(squeeze(fold_TEST(:,:,dimension(j))),  squeeze(fold_TRAIN(:,:,j)) ,fold_TRAIN_GROUP, cfg.classifier_type, observation_weights); % (CVO.training(fi)) for weighted classification and folds

            if strcmp(cfg.classifier_type{1},'toi_classifiers')==1, if strcmp(cfg.classifier_type{2},'trn')==1, txtAuc{fi,j} = hh; continue; end, end % because this classifier needs the time info so we use the variables differently to just get the classifier model at every timept
            if strcmp(cfg.perf_measure,'auc')==1
                % est is now trltime
                forAuc = reshape( hh(:,2),dim(1),[] );
                for timPt = 1:dim(3)
                    [~,~,~,txtAuc_tmp{j}.res(fi,timPt,ch),~] = perfcurve(fold_TEST_GROUP, forAuc(:,timPt) ,2); % trn on rows (j) = yaxis
                    % it was txtAuc(j,timPt,fi) but the structs are faster to
                    % handle in memory instead of carrying the whole 3d with every iteration
                end
            elseif strcmp(cfg.perf_measure,'acc')==1
                forAcc = reshape( outclass,dim(1),[] );
                for timPt = 1:dim(3)
                    txtAuc_tmp{j}.res(fi,timPt,ch) = (sum(fold_TEST_GROUP==forAcc(:,timPt))/length(fold_TEST_GROUP));
                end
            elseif strcmp(cfg.perf_measure,'dval')==1
%                 trl_time_dval{j}.res(fold_dval,:,:) = reshape( max(hh,[],2),dim(1),[] ); % trntime_fold_trl_time
                trl_time_dval{j}.res(fold_dval,:,:) = reshape( outclass,dim(1),[] ); % to return the predicted labels if we want to look at the correct and incorrect trials like in REM_classification_theta
                % every trn timept will give you a whole testing set of fidelity values of trl_time
            end
            %             progressbar.increment();
            %             if strcmp(cfg.classifier_type{1},'toi_classifiers')==1, if strcmp(cfg.classifier_type{2},'tst')==1, break; end, end % because we don't loop on trn time we already got the models at every timept so we just evaluate on all tst pts once
        end
        %         delete(progressbar);
    end

end

if exist('trl_time_dval')==1, txtAuc = cell2mat((cellfun(@(x) (x.res),trl_time_dval,'Un',0 )) ); return;  end % we are returning dval and have trials returned not only time

if strcmp(cfg.classifier_type{1},'toi_classifiers')==1, if strcmp(cfg.classifier_type{2},'trn')==1, return; end, end
% reshaping the result in matrix of classification performance
txtAuc = cell2mat((cellfun(@(x) (mean(x.res,1)),txtAuc_tmp,'Un',0 ))' );
% if we are classifying on every channel they will be occupying the third dimension

end









function [outclass, err, posterior] = low_lvl_classify(TEST,  TRAIN , TRAIN_GROUP, classifier, observation_weights)
% takes the train and test as 2d (trialsxfeatures) and the labels and does classification
% according to the classifier
classifier_name = classifier{1};
rng(1); % For reproducibility
if length(classifier)==3, if strcmp(classifier{2},'domain_align'), [outclass, err, posterior] = lv_txt_domain_adaptation(TRAIN,TEST,TRAIN_GROUP,classifier{3},classifier); return; end,  end% domain alignment of the datasets
% to make sure that the trn boundary should work on test and it should be
% here because every trn time pt has a different classifier and a
% different domain .. classifier{3} has the original dimensions trl_ch_time
switch classifier_name
    case 'lv_lda'
        % FDA matthias' paper
        cfg=[]; cfg.TRAIN=TRAIN; cfg.TEST=TEST; cfg.TRAIN_GROUP=TRAIN_GROUP; cfg.phase='trn'; cfg.mdl=[];
        cfg.tst_noise_for_within=1; % for covariance estimation based on within class covariance/noise
        % of test dataset (may help when noise of trn andd tst are different)
        [mdl] = lv_lda(cfg);
        cfg.mdl=mdl; cfg.phase='predict'; [mdl] = lv_lda(cfg);
        outclass = mdl.outclass; posterior = mdl.posterior; err=1;
    case 'lda'
        %[outclass, err, posterior] = classify(TEST,  TRAIN , TRAIN_GROUP);
        Mdl = fitcdiscr(TRAIN,TRAIN_GROUP,'Weights',observation_weights);
        [outclass,posterior] = predict(Mdl,TEST); err=1;
    case 'svm'
        Mdl = fitcsvm(TRAIN,TRAIN_GROUP, 'KernelFunction' ,'rbf');
        Mdl = fitSVMPosterior(Mdl);  % to fit a sigmoid and get the posterior for the classes o.w. it will be like regression .. you have to make sure.
        [outclass,posterior] = predict(Mdl,TEST); err=1;
    case 'random_forest'
        t = templateTree('NumVariablesToSample','all',...
            'PredictorSelection','interaction-curvature','Surrogate','on');
        % normal training
        Mdl = fitcensemble(TRAIN,TRAIN_GROUP,'Method','Bag','NumLearningCycles',5, ...
            'Learners',t); err = 1;
        [outclass,posterior] = predict(Mdl,TEST);

    case 'naive_bayes'
        [outclass, err, posterior] = my_naive_classifier(TEST,  TRAIN , TRAIN_GROUP);

    case 'non_parametric_kernel_distribution'
        cfg=[]; cfg.TRAIN=TRAIN; cfg.TEST=TEST; cfg.TRAIN_GROUP=TRAIN_GROUP;
        [outclass, posterior] = lv_kernel_distribution(cfg);
        err=1;
    case 'py_classifiers'
        name = classifier{2};
        activation = classifier{3}; % activation .. ex: 'relu'
        %         classifier = py.importlib.import_module(classifier_name); %
        %         doesn't work in matlab2021
        eval(['prediction = classifier.' name '(TEST,  TRAIN , TRAIN_GROUP, activation);']); % name ex: NN
    case 'Riemannian'
        riem_method = classifier{2};
        [outclass, posterior] = lv_Riemannian_classify(TEST,  TRAIN , TRAIN_GROUP, riem_method);
        err=1;
    case 'toi_classifiers'
        if strcmp(classifier{2},'trn')==1 % classifier{2} is phase trn or tst , classifier{3} is the sbj
            % this is training of many classifiers at every timept will put
            % the mdl in the posterior
            Mdl = fitcdiscr(TRAIN,TRAIN_GROUP);
            posterior = Mdl; err=1; outclass=1;
        end
        if strcmp(classifier{2},'tst')==1
            toi_models = classifier{3}; % 1_models(model for every time point in the training set) with additional field in the model called lv_weights which is the weight for every classifier
            w = toi_models.lv_weights; err=1;
            for i=1:size(toi_models.mdl,2)
                [~,posteriors(:,:,i)] = predict(toi_models.mdl{1,i}, TEST);  % fine for lda and some classifiers for svm will need fitPosterior first but lda is fine
                % posteriors: trls_classes_timeptModel
                posteriors(:,:,i) = posteriors(:,:,i) .* w(1,i);
            end
            posterior = sum(posteriors,3);
            [~,outclass] = max(posterior,[],2);
        end
    case 'fully_connected' % network of only a fully connected layer
        TRAIN_GROUP=categorical(TRAIN_GROUP);
        layers = [  % this is a linear preceptron we may add activation here for non-linear classification, examples: leakyReluLayer and sigmoidLayer after featureInputLayer
            featureInputLayer(size(TRAIN,2))
            sigmoidLayer
            fullyConnectedLayer(length(unique(TRAIN_GROUP)))
            softmaxLayer
            classificationLayer];
        options = trainingOptions('adam','InitialLearnRate',0.01, ...
            'MaxEpochs',5,'Verbose',false,'Plots','training-progress'); % ,'Plots','training-progress'
        net = trainNetwork(TRAIN,TRAIN_GROUP',layers,options);
        % NN test
        [YPred] = classify(net,TEST);

        p = str2num(char(YPred));%coverting the result to vector
        outclass = p(:);  err=0;
        % posterior
        layer_name = net.Layers(5, 1).Name;
        posterior = (activations(net, TEST, layer_name,'OutputAs','columns' ) )';

    case {'RNN','TCN','RNN_biLSTM','1d_cnn'} % for deep nets  like RNNs and TCNs
        if strcmp(classifier_name,'RNN')==1
            %         TRAIN2{1} = cell2mat(TRAIN); % to aggregate trials in one long sequence
            %         TRAIN_GROUP2{1}  = categorical(cell2mat(TRAIN_GROUP));
            % nan values are considered missing (jittering) so we remove them
            onePoint = 1;
            if onePoint==1
                % if we have one label and not axtime classification
                for i=1:length(TRAIN_GROUP), ss(1,i) = TRAIN_GROUP{1,i}(1); end
                outcome = 'last';
                TRAIN_GROUP = ss';
            else
                outcome = 'sequence';
                TRAIN_GROUP = cellfun(@(x,y) y(1:length(x)),TRAIN,TRAIN_GROUP, 'Un',0);
            end

            TRAIN = cellfun(@(x) x(:,~isnan(mean(x,1))),TRAIN, 'Un',0);
            TEST = cellfun(@(x) x(:,~isnan(mean(x,1))),TEST, 'Un',0);
            
            % RNN train
            inputSize = size(TRAIN{1},1);
            numHiddenUnits = 50;
            numClasses = 2;%length(unique(cell2mat(cellfun(@(x) (str2double(categories(x(1)))),TRAIN_GROUP,'Un',0))));
            layers = [ ...
                sequenceInputLayer(inputSize)
                lstmLayer(numHiddenUnits,'OutputMode',outcome) %'StateActivationFunction','tanh','GateActivationFunction','sigmoid'
                % bilstmLayer can be used
                fullyConnectedLayer(numClasses)
                softmaxLayer
                classificationLayer];
            options = trainingOptions('adam', ...
                'MaxEpochs',500, ... % was 60 
                'Verbose',0, 'Plots','training-progress'); %...'Plots','training-progress');
            net = trainNetwork(TRAIN,TRAIN_GROUP,layers,options);  % XTrain: guess here it was sbjs for us it will be trials

            %             net.Layers(2, 1).StateActivationFunction='purelin';

            % RNN test
            [YPred] = classify(net,TEST);
            classes_labels = [categorical(cellstr('c1')) ;categorical(cellstr('c2')); ...
                categorical(cellstr('c3')); categorical(cellstr('c4'))];
            if onePoint~=1
                for i=1:size(YPred,1) %trials
                    for j=1:size(YPred{i,1},2) %time
                        temp(i,j) = find(ismember(classes_labels,YPred{i,1}(j)));
                    end
                end
            else
                for i=1:size(YPred,1) %trials
                        temp(i,1) = find(ismember(classes_labels,YPred(i,1)));
                end
                temp = repmat(temp,1,size(TRAIN{1, 1},2));
            end
            %p = cell2mat(cellfun(@(x) str2num(char((x)))',YPred,'Un',0));%coverting the result to vector
            outclass = temp(:); posterior=0; err=0;


        elseif strcmp(classifier_name,'1d_cnn')==1
            % 1D convolution, takes sequence and classficiation o/p is one pt. for sequence
            % to sequence use TCN
            % convolves with the time dimension and then pooling so the
            % time dimension is reduced then it takes the mean of the
            % remaining pts in the time dim. then classifies trial_ch point
            TRAIN = cellfun(@(x) x(:,~isnan(mean(x,1))),TRAIN, 'Un',0);
            TEST = cellfun(@(x) x(:,~isnan(mean(x,1))),TEST, 'Un',0);
            use_another_classifier = {'lda'};
            % 1d_cnn train
            inputSize = size(TRAIN{1},1); 
            numClasses = 4;%length(unique(cell2mat(cellfun(@(x) (str2double(categories(x(1)))),TRAIN_GROUP,'Un',0))));
            filterSize = 20; % was 3 .. 20 for 100ms as each pt is 5ms
            numFilters = 32;

            layers = [ ...
                sequenceInputLayer(inputSize)
                convolution1dLayer(filterSize,numFilters,Padding="causal")
                reluLayer
                layerNormalizationLayer
                convolution1dLayer(filterSize,2*numFilters,Padding="causal")
                reluLayer
                layerNormalizationLayer
                globalAveragePooling1dLayer
                fullyConnectedLayer(numClasses)
                softmaxLayer
                classificationLayer];

            options = trainingOptions('adam', ...
                'MaxEpochs',60, ... 
                'Verbose',0, ... %'Plots','training-progress', ...
                SequencePaddingDirection="left");
            net = trainNetwork(TRAIN,TRAIN_GROUP,layers,options);  % XTrain: guess here it was sbjs for us it will be trials
            
            if ~isempty(use_another_classifier)
                layer_name = net.Layers(8, 1).Name;
                features_TEST = (activations(net, TEST, layer_name,'OutputAs','columns' ) )';
                features_TRAIN = (activations(net, TRAIN, layer_name,'OutputAs','columns' ) )';
                [outclass, err, posterior] = low_lvl_classify(features_TEST,  features_TRAIN , TRAIN_GROUP, use_another_classifier);
                outclass = repmat(str2num(char(outclass)), 1,size(TEST{1,1},2));
            else
                % 1d_cnn test
                [YPred] = classify(net,TEST);
                temp = repmat(str2num(char(YPred)), 1,size(TEST{1,1},2));
                outclass = temp(:); posterior=0; err=0;
            end

        elseif strcmp(classifier_name,'RNN_biLSTM')==1 % different bidirectional network and with batches
            % should use the whole batch  if you can fit it o.w. use
            % batches like 64 and 128
            % RNN train
            inputSize = size(TRAIN{1},1);
            numHiddenUnits = 20;
            numClasses = 2;%length(unique(cell2mat(cellfun(@(x) (str2double(categories(x(1)))),TRAIN_GROUP,'Un',0))));
            miniBatchSize = 64;

            TRAIN = TRAIN(:); TEST = TEST(:); %temp=[];
            %             for i=1:length(TRAIN_GROUP), temp(i,1) = (str2num(char(TRAIN_GROUP{1,i}(1,1)))); end
            %             TRAIN_GROUP = categorical(temp);

            layers = [ ...
                sequenceInputLayer(inputSize)
                bilstmLayer(numHiddenUnits,'OutputMode','sequence')
                fullyConnectedLayer(numClasses)
                softmaxLayer
                classificationLayer];

            maxEpochs = 100;

            options = trainingOptions('adam', ...
                'ExecutionEnvironment','auto', ...
                'GradientThreshold',2, ...
                'MaxEpochs',maxEpochs, ...
                'MiniBatchSize',miniBatchSize, ...
                'Verbose',0, ...
                'Plots','training-progress');

            net = trainNetwork(TRAIN,TRAIN_GROUP,layers,options);
            % RNN test
            YPred = classify(net,TEST, 'MiniBatchSize',miniBatchSize );
            %             p=[];
            %             for i=1:length(YPred), p(i,1) = (str2num(char(YPred(i,1)))); end %coverting the result to vector
            p = cell2mat(cellfun(@(x) str2num(char((x)))',YPred,'Un',0));%coverting the result to vector
            outclass = p(:); posterior=0; err=0;

        elseif strcmp(classifier_name,'TCN')==1 % temporal convolutional network (TCN) routine
            %             [outclass, err, posterior] = lv_TCN_classification(TEST,  TRAIN , TRAIN_GROUP);
            % data is in variable s then we use the script from matlab to build the
            % temporal convolutional network by stacking many CNNs on each other ..
            %creating new trials that have both classes .. the code needs
            %this step to be able to do onehot encoding later on
            lbls=cell2mat(cellfun(@(x) (cat2dbl(x(1))),TRAIN_GROUP,'Un',0));
            idx1 = find(lbls==1); idx2 = find(lbls==2);
            for i=1:min([length(idx1) length(idx2)])
                s.XTrain{i,1} = [single(TRAIN{1,idx1(i)}) single(TRAIN{1,idx2(i)})];
                s.YTrain{i,1} = categorical([repmat(lbls(idx1(i)),1,size(TRAIN{1,idx1(i)},2)) repmat(lbls(idx2(i)),1,size(TRAIN{1,idx2(i)},2))]);
            end
            %             s.XTrain = TRAIN(:); s.XTest = TEST(:);
            %             s.YTrain = TRAIN_GROUP';
            % the loop to min() will not consider the rest of trials because we want to
            % balance the classes so we remove the trials of the class that
            % is repeated more
            s.XTest = TEST';

            %             s.XTrain{1} = single(cell2mat(TRAIN));
            %             s.YTrain{1}=cell2mat(cellfun(@(x) (cat2dbl(x)),TRAIN_GROUP,'Un',0));
            seq2seq_TCN; % follow the procedure of TCN training and testing

        end
    case 'CNN'
        % reshaping back from a vector in the second dimension to the 2d topo at the 2nd dim that it should be
        % the dimesions of topos images are fixed to be: 64*64
        TRAIN = reshape(TRAIN, size(TRAIN,1),sqrt(size(TRAIN,2)),sqrt(size(TRAIN,2)) );
        temp(1,:,:,:) = TRAIN;
        % numeric array of images should have this order of dimensions (hight width no.channels(color) no.images)
        TRAIN = permute(temp,[3 4 1 2]); clear temp;

        %         TRAIN(:,:,:,TRAIN_GROUP==1) = 0; % forcing classes to verify that the network is working
        %         TRAIN(:,:,:,TRAIN_GROUP==2) = 1;


        TRAIN_GROUP=categorical(TRAIN_GROUP);



        layers = [
            imageInputLayer([10 10 1])

            convolution2dLayer(3,8,'Padding','same') % filter_size = 3x3 .. no.filters = 8
            batchNormalizationLayer
            reluLayer

            maxPooling2dLayer(2,'Stride',2) % downsampling 2x2 maximum inside this mask so it's downsampling ..
            % stride is the jumps so we skip a point we don't slide at every pt

            convolution2dLayer(3,16,'Padding','same')
            batchNormalizationLayer
            reluLayer

            maxPooling2dLayer(2,'Stride',2)

            convolution2dLayer(3,32,'Padding','same')
            batchNormalizationLayer
            reluLayer


            fullyConnectedLayer(2)
            softmaxLayer
            classificationLayer];

        options = trainingOptions('adam', ...
            'InitialLearnRate',0.01, ...
            'MaxEpochs',100, ...
            'Verbose',false, ...
            'Plots','training-progress');

        net = trainNetwork(TRAIN,TRAIN_GROUP',layers,options);
        % RNN test
        TEST = reshape(TEST, size(TEST,1),sqrt(size(TEST,2)),sqrt(size(TEST,2)) );
        temp(1,:,:,:) = TEST;
        TEST = permute(temp,[3 4 1 2]); clear temp;
        [YPred] = classify(net,TEST);

        p = str2num(char(YPred));%coverting the result to vector
        outclass = p(:);  err=0;
        % posterior
        layer_name = net.Layers(14, 1).Name; % because 14 is the softmax layer
        posterior = (activations(net, TEST, layer_name,'OutputAs','columns' ) )';
    otherwise
        fprintf(['\n No classifier chosen will assume LDA. \n']);
        [outclass, err, posterior] = classify(TEST,  TRAIN , TRAIN_GROUP);
end


end


function converted=cat2dbl(catdata)
catnames = categories(catdata);
converted = str2double(catnames(catdata)).';
converted=converted';
end



function  [outclass, err, posterior] = my_naive_classifier(TEST,  TRAIN , TRAIN_GROUP)
% Naive bayes classifier
% models the features data into normal distribution and get the posterior
% according to naive theroem .. it's naive because it assumes that the
% features are independent so we get the probabilities of the features
% given the class and multiply then with each other so by that we assume
% they are independent because if they were it wouldn't be just multiplication ..
% p(y|x) = p(x|y) . p(y) / p(x)
% p(x) is constant so it won't change the result

% classes mus and sigmas
mus = [mean(TRAIN(TRAIN_GROUP==1,:),1) ; mean(TRAIN(TRAIN_GROUP==2,:),1)]; % [mu1;mu2]
sigmas = [std(TRAIN(TRAIN_GROUP==1,:),[],1) ; std(TRAIN(TRAIN_GROUP==2,:),[],1)]; % [sigma1;sigma2]
priors = [sum(TRAIN_GROUP==1);sum(TRAIN_GROUP==2)]./size(TRAIN,1); % [priorC1;proprC2]

% posterior for the testing set
for i=1:size(mus,2)
    test_pxy(:,i,1) = normpdf(TEST(:,i),mus(1,i),sigmas(1,i)); % trls_feat_likelihood
    test_pxy(:,i,2) = normpdf(TEST(:,i),mus(2,i),sigmas(2,i));
end


likelihood = [prod([squeeze(test_pxy(:,:,1)) repmat(priors(1), size(test_pxy,1),1) ], 2) ...
    prod([squeeze(test_pxy(:,:,2)) repmat(priors(2), size(test_pxy,1),1) ], 2) ];

posterior = likelihood ./ repmat(sum(likelihood,2),1,size(likelihood,2)); % actual certainty

% % making posterior actual certainty that sums to 1
% post = sum(posterior,2);
% % basically looking at how much they need to the sum to be 1 and then each peobability half of the needed amount
% posterior( post<1,: ) = posterior( post<1,: ) + (1- post(post < 1))./2;
% % if more than 1 then reduce
% posterior( post>1,: ) = posterior( post>1,: ) - abs((1- post(post > 1))./2);


err = 1; % dummy as we don't use it anyway
outclass = (posterior(:,2)>posterior(:,1) )+1; % class 1 or 2

end



function [outclass, posterior] = lv_Riemannian_classify(TEST,  TRAIN , TRAIN_GROUP, method)
% takes TRAIN and TEST they are cov matrices of every trial the cov. is
% reshaped as vector in the second dimension and in here we reshape them
% again to cov and do the classification
% Riemannian geometry classifiers contains three methods:
% (1) distance.. for classifying based on riemannian distance from riemannian mean on
% the manifold (will map the distance to be like posterior)
% (2) tangent .. to classify with actual classifier after projecting the cov.
% data points from the manifold to the tangent space they become
% linear and can use any linear classifier, the tangent space is
% calculated using the mean cov which is calculated with the riem.
% mean.
% (3) pairing .. to classify based on the riem. distance on the manifold
% without calculating the riem mean so that the distances between all
% pair of trials are calculated and then fed to a classifier
% it's fine if zscoring was performed because zscoring is done on
% the cov reshaped in vector so still the same point in the cov
% matrix centered among all trials which is good
ch = sqrt(size(TRAIN,2)); % original channels count because this is cov so chxch so we get sqrt
COVtrain = nan(ch,ch,size(TRAIN,1));
COVtest = nan(ch,ch,size(TEST,1));
for i=1:size(TRAIN,1), COVtrain(:,:,i) = reshape(TRAIN(i,:),ch,ch); end
for i=1:size(TEST,1), COVtest(:,:,i) = reshape(TEST(i,:),ch,ch); end

metric_mean = 'riemann'; % a lot others available: {'euclid','logeuclid','riemann','ld'}
metric_dist = 'riemann'; % a lot others available: {'euclid','logeuclid','riemann','ld'}

Ytrain  = TRAIN_GROUP;

switch method
    case 'distance'
        [outclass,distances] = mdm(COVtest,COVtrain,Ytrain,metric_mean,metric_dist);
        posterior = 1-(distances ./ repmat(sum(distances,2),1,size(distances,2))); % softmax the distance to posterior
    case 'tangent'
        classifier = {'lin_svm'}; % another classifier that operates in the tangent space
        [outclass,posterior] = tslda_matlab(COVtest,COVtrain,Ytrain,metric_mean,metric_dist, classifier);
    case 'pairing'
        classifier = {'lin_svm'}; % another classifier that operates on distances
        [outclass,posterior] = lv_mdm(COVtest,COVtrain,Ytrain,metric_mean,metric_dist, classifier);
end

end




function [outclass, err, posterior] = lv_txt_domain_adaptation(TRAIN,TEST,TRAIN_GROUP,dims,classifier)
% instead of applying the model once here we do it at every timr pt from
% trn and tst bebcause we align the domain so trn and tst will change (not always)
% should take the trialtime_ch data and puts it back into 3d and then
% align every time pt (trialxch) with the TRAIN trialxch then puts it
% back in trialtime_ch for classification ... it aligns both trn and tst so
% it will produce a different trn and tst set at every time pt and then
% calls the lvl classification
% the source is the test set and it aligns to the target which is the training set
[outclass,err,posterior]=deal([]);
source = reshape( TEST,dims(1),dims(3),dims(2));
source = permute(source,[1 3 2]); classifier(2:end)=[];
apply_once=1;

if apply_once==0
    % progressbar = ParforProgressbar(size(source,3),'title', 'Domain align at every time pt.');
    for i=1:size(source,3)
        cfg=[]; transformed=[];
        cfg.trn.trial=squeeze(source(:,:,i)); cfg.tst.trial=TRAIN;
        cfg.method = 'PCA_correlation'; %'domain_align_toolbox';
        cfg.pool = 1;
        if strcmp(cfg.method,'domain_align_toolbox') % because the formatting is different
            % ssa: reduces the features to one stationary feature that maximises the correlation between features
            % may change this to another method inside the toolbox but don't use a
            % method that uses the labels of the test set because we shouldn't use the
            % labels of the test set.. that's why the test target value is set to zero
            % which is the third paramter in ftTrans_ .. change the method to somthing else from the toolbox
            % you may also call a function from lv_component_analysis and
            % transform the tst set by training on a transformation on trn (i.e., PCA)
            ftAll = [cfg.tst.trial ; cfg.trn.trial];
            maSrc = logical([ones(size(cfg.tst.trial,1),1);zeros(size(cfg.trn.trial,1),1)]);
            maLabeled = ones(size(ftAll,1),1);

            [ftAllNew,~] = ftTrans_ssa(ftAll,maSrc,0,maLabeled);

            transformed.target = ftAllNew(1:size(cfg.tst.trial,1),:);
            transformed.source = ftAllNew(size(cfg.tst.trial,1)+1:end,:);
        else
            transformed = lv_feature_extractor(cfg);

            %             cfg=[]; % to add another level of domain adaptation
            %             cfg.trn.trial=transformed.source; cfg.tst.trial=transformed.target;
            %             cfg.method = 'PCA_correlation'; cfg.pool = 0;
            %             transformed = lv_feature_extractor(cfg);

        end

        [outclass_temp, err_temp, posterior_temp] = low_lvl_classify(transformed.source,  transformed.target,...
            TRAIN_GROUP, classifier);
        outclass = [outclass ;outclass_temp]; err = [err ;err_temp]; posterior = [posterior ;posterior_temp];
        %     progressbar.increment();
    end
    % delete(progressbar);
else
    % incase the transformation requires the trn set at specific time pt and this
    %can transform the whole testing set at once like kernel_PCA
    cfg=[]; transformed=[];
    cfg.trn.trial=TEST; cfg.tst.trial=TRAIN;
    cfg.method = 'PCA_correlation'; %'domain_align_toolbox';
    cfg.pool = 0; cfg.labels = TRAIN_GROUP;
    if strcmp(cfg.method,'domain_align_toolbox')
        ftAll = [cfg.tst.trial ; cfg.trn.trial];
        maSrc = logical([ones(size(cfg.tst.trial,1),1);zeros(size(cfg.trn.trial,1),1)]);
        maLabeled = ones(size(ftAll,1),1);

        [ftAllNew,~] = ftTrans_ssa(ftAll,maSrc,0,maLabeled);

        transformed.target = ftAllNew(1:size(cfg.tst.trial,1),:);
        transformed.source = ftAllNew(size(cfg.tst.trial,1)+1:end,:);
    else
        % %         % for CNN
        % %         cfg=[]; cfg.trn.trial=TRAIN; cfg.tst.trial=TEST;
        % %         cfg.trn.trialinfo(:,1) = TRAIN_GROUP;
        % %         transformed_temp = lv_EEG_to_CNN(cfg);
        % %         transformed.source = transformed_temp.tst{1, 1}.trial; transformed.target = transformed_temp.trn.trial;

        transformed = lv_feature_extractor(cfg);

        %         cfg=[]; % to add another level of domain adaptation
        %         cfg.trn.trial=transformed.source; cfg.tst.trial=transformed.target;
        %         cfg.method = 'PCA_correlation'; cfg.pool = 0;
        %         transformed = lv_feature_extractor(cfg);
    end
    [outclass, err, posterior] = low_lvl_classify(transformed.source,  transformed.target,...
        TRAIN_GROUP, classifier);

end


end


function [outclass, posterior] = lv_kernel_distribution(cfg)
% kernel non parametric distribution
% uses the kernel technique which assumes that every data point is a
% center of a normal distribution and then we use a std. h or bandwidth
% to control the smoothing of every distribution (the higher the more smoothing) and them we sum
% all the distributions and normalise into one distribution which is
% now the pdf and can be used to project a new data point .. in matlab
% they have a function ready for this and the great thing about it that
% it enables changing the std. for different features .. but still
% kernel method is controlled by the std. and assumes that the density
% of data pts at different values are with the same std.
h = std(cfg.TRAIN,[],1);
d = size(cfg.TRAIN,2); n = size(cfg.TRAIN,1);
bw = h * (4./(n*(d+2))).^(1./(d+4)); %  Silverman, B. W. Density Estimation for Statistics and Data Analysis. Chapman & Hall/CRC, 1986.

id1 = cfg.TRAIN_GROUP==1; id2 = cfg.TRAIN_GROUP==2;
likelihood1 = mvksdensity(cfg.TRAIN(id1,:),cfg.TEST,'Bandwidth',bw);
likelihood2 = mvksdensity(cfg.TRAIN(id2,:),cfg.TEST,'Bandwidth',bw);
% softmax the likelihood to get the posterior
posterior = [likelihood1./(likelihood1+likelihood2)  likelihood2./(likelihood1+likelihood2)];
outclass = 1*(posterior(:,1)>posterior(:,2))+2*(posterior(:,1)<posterior(:,2));
end


% function [new_TRAIN,new_TEST] = lv_call_align(TRAIN,TEST,dims)
% % should take the trialtime_ch data and puts it back into 3d and then
% % align every time pt (trialxch) with the TRAIN trialxch then puts it
% % back in trialtime_ch for classification
% source = reshape( TEST,dims(1),dims(3),dims(2)); cfg=[]; new_source=[];
% source = permute(source,[1 3 2]);
% for i=1:size(source,3)
%     cfg.trn.trial=squeeze(source(:,:,i)); cfg.tst.trial=TRAIN;
%     cfg.method = 'domain_align_toolbox'; if strcmp(cfg.method,'domain_align_toolbox'), cfg.trn.trial = [cfg.tst.trial ; cfg.trn.trial]; end
%     transformed = lv_feature_extractor(cfg);
%     if strcmp(cfg.method,'domain_align_toolbox') % will only work for methods where we don't transform the target TRAIN..
%         %transformed.target=transformed.target(1:size(cfg.tst.trial,1),:);
%         transformed.source=transformed.source(size(cfg.tst.trial,1)+1:end,:);
%     end
%     new_source = [new_source ; transformed.source];
% end
% new_target = transformed.target; % this happens once because the target trn doesn't change in
% % SA, CORAL, and EA,in EA it changes but not relative to tst so we do this
% % change once in here not with every tst timept.
%
% new_TEST = new_source; new_TRAIN = new_target;
% end
% % function [TRAIN,TEST] = lv_domain_adaptation(TRAIN, TEST)
% % source = reshape( TEST,dims(1),dims(3),dims(2)); cfg=[]; new_source=[];
% % source = permute(source,[1 3 2]);
% % for i=1:size(source,3)
% %     cfg.trn.trial=squeeze(source(:,:,i)); cfg.tst.trial=TRAIN;
% %     ftAll = [TRAIN ; squeeze(source(:,:,i))];
% %
% %     maSrc = logical([ones(size(TRAIN,1),1);zeros(size(squeeze(source(:,:,i)),1),1)]);
% %     maLabeled = ones(size(ftAll,1),1);
% %     % trying SSA domain adaptation ..
% %     % reducing the features to one stationary feature that maximises the correlation between features
% %     % may change this to another method inside the toolbox but don't use a
% %     % method that uses the labels of the test set because we shouldn't use the
% %     % labels of the test set.. that's why the test target value is set to zero
% %     % which is the third paramter in ftTrans_
% %     [ftAllNew,~] = ftTrans_ssa(ftAll,maSrc,0,maLabeled);
% %     % another example: [ftAllNew,transMdl] =
% %     % ftTrans_pca(ftAll,maSrc,target,maLabeled,param); ftAllNew = ftAllNew(:,1:3); % 3 PCs
% %
% %     TRAIN = ftAllNew(1:size(TRAIN,1),:); new_source = [new_source ; ftAllNew(size(TRAIN,1)+1:end,:)];
% %
% % end
% % new_target = TRAIN; % this happens once because the target trn doesn't change in
% % % SA, CORAL, and EA,in EA it changes but not relative to tst so we do this
% % % change once in here not with every tst timept.
% % TEST = new_source; TRAIN = new_target;
% %
% %
% % end

% optimized non parallel
%
% for j=1:size(fold_TRAIN,3)
%         % always take the first slice and shrink the time to get faster as we reach more iterations (fold_TRAIN size decreases)
%         tic
%         [~,~,hh] = low_lvl_classify(squeeze(fold_TEST(:,:,1)),  squeeze(fold_TRAIN(:,:,1)) ,fold_TRAIN_GROUP, cfg.classifier_type);
%         toc
%         tic
%         % est is now trltime
%         forAuc = reshape( hh(:,2),dim(1),[] );
%         for timPt = 1:dim(3)
%             [~,~,~,txtAuc_tmp{j}.res(fi,timPt),~] = perfcurve(fold_TEST_GROUP, forAuc(:,timPt) ,2); % trn on rows (j) = yaxis
%             % it was txtAuc(j,timPt,fi) but the structs are faster to
%             % handle in memory instead of carrying the whole 3d with every iteration
%         end
%         toc
%         %progressbar.increment();
%         fold_TRAIN(:,:,1)=[];
%         if size(fold_TEST,3)~=1, fold_TEST(:,:,1)=[]; end % if the data is axtime we will delete the slice on top and
%         % proceed on time with every new iteration
%     end




