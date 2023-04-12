% BCI competition IV dataset 1
% Benjamin Blankertz, Guido Dornhege, Matthias Krauledat, Klaus-Robert Müller, and Gabriel Curio. The non-invasive
% Berlin Brain-Computer Interface: Fast acquisition of effective performance in untrained subjects. NeuroImage,
% 37(2):539-550, 2007.

clear all
%load data and format
ppnt_names = {'b','c','d','e','g'}; % removed a and f because not left vs right hand
raw_path = '.\data\BCICIV_calib_ds1';
TF_temp=[];

% method
power_classification_axtime = 0;
csp = 0;
riemannian = 0;
RNN = 1;
%%
for nn=1:length(ppnt_names)
    load([raw_path ppnt_names{nn} '.mat']);
    dat.label = nfo.clab';
    dat.trial{1,1} = 0.1*double(cnt');
    dat.sampleinfo = [1 size(dat.trial{1,1},2)];
    dat.fsample = nfo.fs;
    dat.time{1,1} = 0:1/dat.fsample:(dat.sampleinfo(2)-1)/dat.fsample;

    %% filtering the continuous data
    cfg=[]; %the 100hz version was already band pass filtered and band stop
    %filtered
    %cfg.bpfilter  = 'yes'; % because the data was already filtered 0.05 and 200 Hz
    %cfg.bpfreq    = [0.1 50]; % for new data sets consider FIR because it's better for offline analyses and not data driven as butter and IIR !
    %cfg.bsfilter  = 'yes'; % band-stop method
    %cfg.bsfreq    = [48 52];
    cfg.demean      = 'yes';
    dat  = ft_preprocessing(cfg, dat); % Matlab: fieldtrip, EEGlab ... Python: MNE
    cfg=[];
    cfg.channel=lower({'C6','C4','C2','C1','C3','C5', 'CP5', 'CP3', 'CP1', 'CP2', 'CP4', 'CP6'});
    dat = ft_selectdata(cfg, dat);

    %% segmenting into trials
    % because they rested between trials we will take the onset
    % of the cue and segment based on that ... because we extend 4sec. after an onset
    data=[];
    trl_len_sec = 4;
    onsets = mrk.pos'; onsets = repmat(onsets,1,trl_len_sec*dat.fsample+1);
    onsets = onsets + repmat(0:(size(onsets,2)-1),size(onsets,1),1);
    for i=1:size(dat.trial{1,1},1) % looping on channels and putting the trlxtime matrix
        temp= dat.trial{1,1}(i,:);
        data.trial(:,i,:) = temp(onsets);
    end
    data.time = 0:1/dat.fsample:trl_len_sec;
    data.label = dat.label;
    data.sampleinfo = [mrk.pos' (mrk.pos')+(4*dat.fsample)];
    data.trialinfo = mrk.y;
    data.trialinfo(data.trialinfo==-1)=2; % to make classes 1 and 2. it was -1 for class one or 1 for class two so left hand is actually 2 now
    data_trialinfo = data.trialinfo;


    data.trialinfo = data_trialinfo;

    %% visualising the power across time
    if power_classification_axtime==1
        % baseline
        data.baseline=[];
        baseline=3; onsets = mrk.pos'; onsets = onsets-baseline*dat.fsample; % for baseline we go back a bit
        onsets = repmat(onsets,1,baseline*dat.fsample+1);
        onsets = onsets + repmat(0:(baseline*dat.fsample),size(onsets,1),1);
        for i=1:size(dat.trial{1,1},1) % looping on channels and putting the trlxtime matrix
            temp= dat.trial{1,1}(i,:);
            data.baseline(:,i,:) = temp(onsets);
        end
        % to get power across time and see the difference temporally in mu rhythm
        cfg=[]; cfg.data = data; cfg.method = 'power';
        features_data = lv_feature_extractor(cfg);
        data.trial = features_data.trial;
        cfg=[];
        baselinedat = data;  baselinedat.trial = data.baseline; baselinedat.baseline=[];
        baselinedat.time = -baseline:1/dat.fsample:0;
        cfg.data = baselinedat; cfg.method = 'power';
        data.baseline = lv_feature_extractor(cfg);
        baseline_mat(:,:,1) = mean(data.baseline.trial,3);
        baseline_ready = repmat(baseline_mat,1,1,size(data.trial,3));
        data.trial = ((data.trial - baseline_ready)./baseline_ready)*100; % percentage change
        right_hemi_c1(nn,:) = squeeze(median(data.trial(data.trialinfo==2,ismember(data.label,'C4'),:), 1)); % to see power change from baseline axtime .. plotting is done outside the loop to
        % consider all participants
        left_hemi_c1(nn,:) = squeeze(median(data.trial(data.trialinfo==2,ismember(data.label,'C3'),:), 1));

        right_hemi_c2(nn,:) = median(data.trial(data.trialinfo==1,ismember(data.label,'C4'),:), 1);
        left_hemi_c2(nn,:) = median(data.trial(data.trialinfo==1,ismember(data.label,'C3'),:), 1);

        % classification across time
        cfg=[];
        cfg.method = 'axtime'; cfg.classifier_type = {'lda'};
        cfg.do_parallel = 0;
        cfg.trn = data; cfg.trn.trialinfo=data.trialinfo';  cfg.folds=5;
        auc(nn,:) = lv_classify(cfg);
        continue;
    end
    %% classification after CSP
    % with the freq range as a parameter ..
    % CSP is supervised so we need to prevent info. leakage and need to use
    % only the training set to build the csp filters so we will split the
    % data into two sets training and testing ..
    [datatrn,datatst] = deal(data);
    rng(1), idx = randperm(size(data.trial,1));
    idx_trn = idx(1:round(size(data.trial,1)/2)); idx_tst = idx(round(size(data.trial,1)/2)+1:size(data.trial,1));

    datatrn.trial = data.trial(idx_trn,:,:); datatst.trial = data.trial(idx_tst,:,:);
    datatrn.trialinfo = data.trialinfo(idx_trn); datatst.trialinfo = data.trialinfo(idx_tst);
    datatrn_raw = datatrn; datatst_raw = datatst;

    % filtering
    warning('be careful of edges'); % filtering
    cfg_preprocessing                 = [];
    cfg_preprocessing.bpfilter        = 'yes';
    cfg_preprocessing.bpfreq          = [8 13];
    data_bp= ft_preprocessing(cfg_preprocessing, datatrn);
    hold_trn = [log(var(data_bp.trial(:,ismember(data_bp.label,'C3'),:),[],3)) log(var(data_bp.trial(:,ismember(data_bp.label,'C4'),:),[],3))];
    riem_trn = data_bp;

    cfg=[]; % building CSP filters
    cfg.method='csp';
    cfg.step='calculate';
    cfg.data=data_bp;
    cfg.data.trialinfo=datatrn.trialinfo';
    comp = lv_component_analysis(cfg);
    comp.id = [comp.id(1) ;comp.id2(1)]; % choosing csp filters from sorted
    cfg.step='transform'; % applying the chosen filters to transform data
    cfg.comp=comp;
    datatrn.trial = lv_component_analysis(cfg);

    % apply to test
    cfg_preprocessing                 = [];
    cfg_preprocessing.bpfilter        = 'yes';
    cfg_preprocessing.bpfreq          = [8 13];
    data_bp= ft_preprocessing(cfg_preprocessing, datatst);
    hold_tst = [log(var(data_bp.trial(:,ismember(data_bp.label,'C3'),:),[],3)) log(var(data_bp.trial(:,ismember(data_bp.label,'C4'),:),[],3))];
    cfg.data=data_bp;
    riem_tst = data_bp;

    cfg.step='transform';
    cfg.comp=comp;
    datatst.trial = lv_component_analysis(cfg);

    
    if csp==1
        % classification after getting the var over time so the
        % classification would be one point
        % var
        datatrn.trial = log(var(datatrn.trial,[],3));
        datatst.trial = log(var(datatst.trial,[],3));
        datatrn.trial = zscore(datatrn.trial, [], 1); datatst.trial = zscore(datatst.trial, [], 1);
        cfg=[];
        cfg.method = 'axtime'; cfg.classifier_type = {'lda'};
        cfg.do_parallel = 0;
        cfg.trn = datatrn; cfg.trn.trialinfo=datatrn.trialinfo';  cfg.folds=nan;
        cfg.tst = datatst; cfg.tst.trialinfo=datatst.trialinfo';
        auc(nn,:) = lv_classify(cfg);
        continue;
    end


            
    %% classification with Riemannian based geometry on the covariance matrices
    if riemannian==1
        %% extracting covs
        [feat_trn,feat_tst] = deal([]);
        for i=1:size(datatrn.trial,1), feat_trn(i,:,1) = reshape(cov(squeeze(riem_trn.trial(i,:,:))'),1,[]); end
        for i=1:size(datatst.trial,1), feat_tst(i,:,1) = reshape(cov(squeeze(riem_tst.trial(i,:,:))'),1,[]); end
        datatrn.trial = feat_trn; datatst.trial = feat_tst;
        %% Riemannian classification 
        cfg=[];
        cfg.method = 'axtime'; cfg.classifier_type = {'Riemannian','distance'}; % distance tangent pairing
        cfg.do_parallel = 0; cfg.perf_measure='auc';
        cfg.trn = datatrn; cfg.trn.trialinfo=datatrn.trialinfo';  cfg.folds=nan;
        cfg.tst = datatst; cfg.tst.trialinfo=datatst.trialinfo';
        auc(nn,:) = lv_classify(cfg);
        continue;
    end

    %% classification with Recurrent Neural Network
    if RNN==1
        datatrn.trial = datatrn_raw.trial;
        datatst.trial = datatst_raw.trial;

        cfg=[];
        cfg.method = 'axtime'; cfg.method = 'deep'; % distance tangent pairing
        cfg.classifier_type = {'RNN'};
        cfg.do_parallel = 0; cfg.perf_measure='acc';
        cfg.trn = datatrn; cfg.trn.trialinfo=datatrn.trialinfo';  cfg.folds=nan;
        cfg.tst = datatst; cfg.tst.trialinfo=datatst.trialinfo';
        auc(nn,:) = lv_classify(cfg);
        continue;
    end
end

if power_classification_axtime==1
    % visualising and doing the stats. across time to search for
    % significant clusters (with other data we would need more participants in general. e.g., more than 10)
    lv_pretty_errorbar(data.time,left_hemi_c1,left_hemi_c2, 2); figure, % the third parameter is for the stats: 0 for points wise comparison using Wilcoxon, 1 for parametric cluster-based permutation
    % 2 for non parametric cluster-based permutation, 99 for non 
    lv_pretty_errorbar(data.time,right_hemi_c1,right_hemi_c2, 2);
end

auc
mean_auc = mean(auc)

% group lvl TF
TF_struct.trial = TF_temp;
TF_struct.parametric = 1;
lv_tf(TF_struct, 1, 1); %data, do_stats, do_plot


x = randn(10,10);
[V,D] = eig(x)
[U,S,V] = svd(x);



lv_pretty_errorbar(data.time,auc,(auc.*0)+0.5, 2);


lv_pretty_errorbar( auc,(auc.*0)+0.5 );






% hilbert power in mu
figure,
subplot(121), lv_pretty_errorbar(data.time,squeeze(left_hemi_c1),squeeze(left_hemi_c2), 99);
subplot(122), lv_pretty_errorbar(data.time,squeeze(right_hemi_c1),squeeze(right_hemi_c2), 99);
legend([nfo.classes(1) nfo.classes(1) nfo.classes(2) nfo.classes(2)])





%% testing 2d cluster stats
load ('D:\work\result_exp_var_pw')
load ('D:\work\result_adp_var_pw')
vsChance = 1; % vs. chance level (0.5):1  or  the control night: 0

result_exp = result_exp_var_pw;
result_adp = result_exp_var_pw;
times.timeTrn = 0:1/200:1.5;
times.timeTst = 0:1/200:1.1;
% provided as .mat with data
dat=[];
for i=1:size(result_exp,1), cond1(i,:,:)=result_exp{i,1}.perf; cond2(i,:,:)=result_adp{i,1}.perf; end
if vsChance==1, cond2 = (cond2.*0)+0.50; end
id = mod(1: size(cond1,1)+size(cond2,1), 2);
dat.trial(id==0,1,:,:) = cond2; % chance: (cond2.*0)+0.5;
dat.trial(id==1,1,:,:) = cond1;
dat.time=times.timeTst; dat.freq=times.timeTrn; dat.label={'z-stat'};
dat.parametric = 0;
dat.powspctrm=dat.trial;

tic
temp = lv_tf(dat,1,1); %data, do_stats, do_plot
toc





% helping function
function cov_data_pooled=pooled_cov(data, option)
% gets the pooled covariance of different trials (trials in cells)
% riemannian or normal averaging of covs
for i=1:size(data,1), temp{1,i} = squeeze(data(i,:,:)); end
data = temp;
data_centered = cellfun(@(x) (x-repmat(mean(x,2),1,size(x,2)) ), data,'Un',0);
cov_data = cellfun(@(x) ((x*x')./size(x,2)), data_centered,'Un',0);
cov_data_pooled=zeros(size(cov_data{1,1}));

if strcmp(option,'normal')==1
    for i=1:length(data), cov_data_pooled = cov_data_pooled + cov_data{1,i}; end
    cov_data_pooled = cov_data_pooled./length(data);
elseif strcmp(option,'riemann')==1
    covariances = nan( size(data{1,1},1),size(data{1,1},1),length(data) );
    for i=1:length(data), covariances(:,:,i) = cov_data{1,i}; end
    cov_data_pooled = mean_covariances(covariances,'riemann');
end


end