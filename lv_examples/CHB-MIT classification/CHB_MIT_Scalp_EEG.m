% CHB-MIT Scalp EEG Database
% https://physionet.org/content/chbmit/1.0.0/
% we use small subset of the data to do this classification
db_path = '.\data\';
db_name = {'chb01_03' 'chb01_04' 'chb01_15' 'chb01_16' 'chb01_18' 'chb01_21' 'chb01_26'};
startTime = [2996 1467 1732 1015 1720 327 1862];
endTime = [3036 1494 1772 1066 1810 420 1963];
all_features = [];
all_labels = [];

resample = 0; % whether to resample the data to a different sampling rate
cfgResample.resamplefs = 50; % the new sampling rate

% formatting and segmenting data
for nn=1:length(db_name)
    % define trials
    cfg            = [];
    cfg.dataset    = [db_path char(db_name(nn)) '.edf'];
    cfg.continuous = 'yes';
    cfg.channel    = 'all';
    data           = ft_preprocessing(cfg);
    fs=data.hdr.Fs;
    % time is in hrs:minutes:seconds
    % the time in the time vector starts from 0
    % checking that the time in the .txt match the length of the signals in .edf
    % t1 = datetime(2020,1,1,13,43,04,0);
    % t2 = datetime(2020,1,1,14,43,04,0);
    % dt = between(t1,t2)
    if resample == 1
        [data] = ft_resampledata(cfgResample, data);
        fs = cfgResample.resamplefs;
    end

    % cutting into trials
    trl=[];
    offset    = 0;  % the beginning of the trial in case we wanted to include periods before the beginning (baseline)
    trlbegin  = startTime(nn) * fs;
    trlend    = endTime(nn) * fs;
    newtrl    = [trlbegin trlend offset];
    cfg     = []; % getting the trial info
    cfg.trl = newtrl;
    data_def = ft_definetrial(cfg);

    data_class2 = ft_redefinetrial(data_def, data);

    % cutting a random part to represent class1
    temp = data; % taking out the parts with class2
    temp.trial{1,1}(:,trlbegin:trlend) = [];
    rng(1), % for reproducability
    random_start = randi([100 size(temp.trial{1,1},2) - (length(data_class2.time{1,1})+100)],1,1); % take a segment from any part until the end and we cut 100 samples to avoid edges
    dataTOChooseFrom = temp.trial{1,1};
    data_class1 = data_class2;
    data_class1.trial{1,1} = dataTOChooseFrom(:, random_start: random_start+length(data_class2.time{1,1})-1);
    clear temp;

    % looking at the data
    browser_cfg=[];
    browser_cfg.ylim=[-600 600]; % amplitude limit
    browser_cfg.blocksize = 40; % duration for cutting the displayed data, each segment will be 40 seconds.
    browser_cfg.continuous  = 'yes';
    ft_databrowser(browser_cfg,data_class1);
    ft_databrowser(browser_cfg,data_class2);

    % putting the data into trials x channels x time .. each trial will be 2 seconds
    trial_length = 2;
    max_trls = floor(length(data_class1.time{1,1}) / (trial_length*fs));
    trls_idx = 1 : max_trls * (trial_length * fs); % convert the max no. trials into samples to know where to cut the data into equal length trials
    trls_idx = reshape(trls_idx, 2*fs, max_trls);
    trls_idx = trls_idx';

    dataC1=[];
    dataC2=[];

    for ch=1:length(data_class1.label)
        temp = data_class1.trial{1,1}(ch,:);
        dataC1(:,ch,:) = temp(trls_idx);
        temp = data_class2.trial{1,1}(ch,:);
        dataC2(:,ch,:) = temp(trls_idx);
    end


    % getting the fourier spectrum of each trial
    cfg = [];
    cfg.data.trial = dataC1;
    cfg.method = 'spectrum'; cfg.sampling_rate=fs;
    specturms_c1_temp = lv_feature_extractor(cfg);
    specturms_c1 = specturms_c1_temp.trial;

    cfg.data.trial = dataC2;
    specturms_c2_temp = lv_feature_extractor(cfg);
    specturms_c2 = specturms_c2_temp.trial;

    features_c1 = squeeze(mean(specturms_c1,2));
    features_c2 = squeeze(mean(specturms_c2,2));

    % looking at the features
    figure, imagesc(features_c1), caxis([0 100]), title('class1 features')
    figure, imagesc(features_c2), caxis([0 100]), title('class2 features')

    % storing features to aggregate parts
    all_features = [all_features ; features_c1 ; features_c2];
    all_labels = [all_labels ; ones(size(dataC1,1),1) ; 2*ones(size(dataC2,1),1)];

end

% save all_features all_features
% save all_labels all_labels

%%
% load all_features all_features
% load all_labels all_labels


% classification using fully connected layer 
cfg=[]; cfg.method = 'timextime';
cfg.classifier_type = {'fully_connected'}; 
cfg.perf_measure = 'auc';
cfg.trn.trial = all_features; cfg.trn.trialinfo=all_labels; 
cfg.folds=5; cfg.do_parallel=1;
result = lv_classify(cfg); % acc at every trn time
mean(result)


 



