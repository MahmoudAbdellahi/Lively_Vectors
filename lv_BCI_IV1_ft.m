% BCI competition IV dataset 1
% Benjamin Blankertz, Guido Dornhege, Matthias Krauledat, Klaus-Robert Müller, and Gabriel Curio. The non-invasive
% Berlin Brain-Computer Interface: Fast acquisition of effective performance in untrained subjects. NeuroImage,
% 37(2):539-550, 2007.

clear all
%load data and format
ppnt_names = {'b','c','d','e','g'}; % removed a and f because not left vs right hand
raw_path = 'E:\Lively Vectors\BCI_IV\Raw and filtered\BCICIV_calib_ds1';
TF_temp=[];
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
    %     %% this data is filtered but we could say that we need further pre-processing and look at more cleaning steps
    %     data_cleaned = lv_clean_segmented(data, [data.time(1) data.time(end)], nn);
    %
    %     %% manual cleaning
    %     data = data_cleaned;
    %     data.trialinfo = [data.trialinfo' (1:length(data.trialinfo))' ...
    %         data.sampleinfo ones(length(data.trialinfo),1)];
    %
    %
    %     lv_manual_cleaning


    %% ERP analysis on simulated data
% %         signals = lv_drawing_to_signals(100); plot(1:100, signals); % group lvl stats on the erp as if they are coming from different ppnts
% %         signals2 = lv_drawing_to_signals(100); plot(1:100, signals2);
% %         signals = [signals; signals2];
% %     %%
% %         clear temp
% %         load lv_layout.mat
% %         for i=1:length(lv_layout.label), temp(:,i,:) = signals; end
% %         sdata.trial = temp; % adding fields
% %         sdata.trialinfo = [ones(size(signals,1)/2,1);ones(size(signals,1)/2,1)*2];
% %         sdata.time = 0:1/100:1-(1/100);
% %         sdata.layout = lv_layout; sdata.label = lv_layout.label;
% %     
% %         erps_temp = [];
% %         erps = lv_erp(sdata, 0, 0); %data, do_stats, do_plot ... returns 2_ch_time the first dim is cond1 erp then cond2
% %         erps_temp = [erps_temp ; erps.trial]; % erps_temp aggregates all the erps from different sbj
% %         % group lvl ERP
% %         erps.trial = repmat(erps_temp, 30,1);
% %         erps.parametric = 1;
% %         lv_erp(erps, 1, 0); %data, do_stats, do_plot

% %     % check deep with this data
% %     c1(:,1,:) = repmat(signals(1:5,:),100,1,1); c1 = repmat(c1,1,59,1);
% %     c2(:,1,:) = repmat(signals(6:10,:),100,1,1); c2 = repmat(c2,1,59,1);
% %     X_train = permute([c1;c2], [1 3 2]);
% %     save X_train X_train
% %     y_train = [ones(500,1); zeros(500,1) ]; save y_train y_train;
% %     X_test = permute([c1;c2], [1 3 2]);
% %     save X_test X_test
% %     y_test = [ones(500,1); zeros(500,1) ]; save y_test y_test;
    

    
    %% TF analysis
        % baseline
        baseline=3; onsets = mrk.pos'; onsets = onsets-baseline*dat.fsample; % for baseline we go back a bit
        onsets = repmat(onsets,1,baseline*dat.fsample);
        onsets = onsets + repmat(0:(baseline*dat.fsample)-1,size(onsets,1),1);
        for i=1:size(dat.trial{1,1},1) % looping on channels and putting the trlxtime matrix
            temp= dat.trial{1,1}(i,:);
            baseline_dat(:,i,:) = temp(onsets);
        end
        data.trial = cat(3, baseline_dat, data.trial); % aggregating baseline and data
        data.baseline = -baseline:1/dat.fsample:0-1/dat.fsample;
        data.time = [data.baseline data.time];
    
        data.trialinfo = data.trialinfo(:);
        TF_struct = lv_tf(data, 0, 0); %data, do_stats, do_plot .. gets the TF in TF_struct.trial
        TF_temp = [TF_temp ; TF_struct.trial]; % erps_temp aggregates all the erps of different sbj
    
    
        continue

    %% visualising the power across time
    %     % baseline
    %     data.baseline=[];
    %     baseline=3; onsets = mrk.pos'; onsets = onsets-baseline*dat.fsample; % for baseline we go back a bit
    %     onsets = repmat(onsets,1,baseline*dat.fsample+1);
    %     onsets = onsets + repmat(0:(baseline*dat.fsample),size(onsets,1),1);
    %     for i=1:size(dat.trial{1,1},1) % looping on channels and putting the trlxtime matrix
    %         temp= dat.trial{1,1}(i,:);
    %         data.baseline(:,i,:) = temp(onsets);
    %     end
    %     % to get power across time and see the difference temporally in mu rhythm
    %     cfg=[]; cfg.data = data; cfg.method = 'power';
    %     features_data = lv_feature_extractor(cfg);
    %     data.trial = features_data.trial;
    %     cfg=[];
    %     baselinedat = data;  baselinedat.trial = data.baseline; baselinedat.baseline=[];
    %     baselinedat.time = -baseline:1/dat.fsample:0;
    %     cfg.data = baselinedat; cfg.method = 'power';
    %     data.baseline = lv_feature_extractor(cfg);
    %     baseline_mat(:,:,1) = mean(data.baseline.trial,3);
    %     baseline_ready = repmat(baseline_mat,1,1,size(data.trial,3));
    %     data.trial = ((data.trial - baseline_ready)./baseline_ready)*100; % percentage change
    %     right_hemi_c1(nn,:) = squeeze(median(data.trial(data.trialinfo==2,ismember(data.label,'C4'),:), 1)); % to see power change from baseline axtime
    %     left_hemi_c1(nn,:) = squeeze(median(data.trial(data.trialinfo==2,ismember(data.label,'C3'),:), 1));
    %
    %     right_hemi_c2(nn,:) = median(data.trial(data.trialinfo==1,ismember(data.label,'C4'),:), 1);
    %     left_hemi_c2(nn,:) = median(data.trial(data.trialinfo==1,ismember(data.label,'C3'),:), 1);
    %
    %     % classification across time
    %     cfg=[];
    %     cfg.method = 'axtime'; cfg.classifier_type = {'lda'};
    %     cfg.do_parallel = 0;
    %     cfg.trn = data; cfg.trn.trialinfo=data.trialinfo';  cfg.folds=5;
    %     auc(nn,:) = lv_classify(cfg);

    %% feature extraction
    % with the freq range as a parameter ..
    % CSP is supervised so we need to prevent info. leakage and need to use
    % only the training set to build the csp filters so we will split the
    % data into two sets training and testing ..
    [datatrn,datatst] = deal(data);
    rng(1), idx = randperm(size(data.trial,1));
    idx_trn = idx(1:round(size(data.trial,1)/2)); idx_tst = idx(round(size(data.trial,1)/2)+1:size(data.trial,1));

    datatrn.trial = data.trial(idx_trn,:,:); datatst.trial = data.trial(idx_tst,:,:);
    datatrn.trialinfo = data.trialinfo(idx_trn); datatst.trialinfo = data.trialinfo(idx_tst);

    

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


    
 

%         % var or non var
%         % var
%         datatrn.trial = log(var(datatrn.trial,[],3));
%         datatst.trial = log(var(datatst.trial,[],3));
%     
%     
%         % classification at a time
%         datatrn.trial = zscore(datatrn.trial, [], 1); datatst.trial = zscore(datatst.trial, [], 1);
%         cfg=[];
%         cfg.method = 'axtime'; cfg.classifier_type = {'lda'};
%         cfg.do_parallel = 0;
%         cfg.trn = datatrn; cfg.trn.trialinfo=datatrn.trialinfo';  cfg.folds=nan;
%         cfg.tst = datatst; cfg.tst.trialinfo=datatst.trialinfo';
%         auc(nn,:) = lv_classify(cfg);

%     %% extracting covs 
%     [feat_trn,feat_tst] = deal([]);
%     for i=1:size(datatrn.trial,1), feat_trn(i,:,1) = reshape(cov(squeeze(datatrn.trial(i,:,:))'),1,[]); end
%     for i=1:size(datatst.trial,1), feat_tst(i,:,1) = reshape(cov(squeeze(datatst.trial(i,:,:))'),1,[]); end 
%     datatrn.trial = feat_trn; datatst.trial = feat_tst;
%     %% Riemannian classification 
% %     datatrn.trial = zscore(datatrn.trial, [], 1); datatst.trial = zscore(datatst.trial, [], 1);
%     cfg=[];
%     cfg.method = 'axtime'; cfg.classifier_type = {'Riemannian','tangent'}; % distance tangent pairing
%     cfg.do_parallel = 0; cfg.perf_measure='auc';
%     cfg.trn = datatrn; cfg.trn.trialinfo=datatrn.trialinfo';  cfg.folds=nan;
%     cfg.tst = datatst; cfg.tst.trialinfo=datatst.trialinfo';
%     auc(nn,:) = lv_classify(cfg);

    
    %% deep 
    datatrn.trial = datatrn.trial.^2;
    datatst.trial = datatst.trial.^2;
    datatrn.trial = zscore(datatrn.trial, [], 1); datatst.trial = zscore(datatst.trial, [], 1);
    cfg=[];
    cfg.method = 'axtime'; cfg.method = 'deep'; % distance tangent pairing
    cfg.classifier_type = {'RNN'};
    cfg.do_parallel = 0; cfg.perf_measure='auc';
    cfg.trn = datatrn; cfg.trn.trialinfo=datatrn.trialinfo';  cfg.folds=nan;
    cfg.tst = datatst; cfg.tst.trialinfo=datatst.trialinfo';
%     auc(nn,:) = lv_classify(cfg);


% preparing data to send it to keras
    X_train = permute(datatrn.trial,[1 3 2]);
    y_train = datatrn.trialinfo';
    X_test = permute(datatst.trial,[1 3 2]);
    y_test = datatst.trialinfo';
    y_test(y_test==2)=0; y_train(y_train==2)=0;
    save X_train X_train, save y_train y_train 
    save X_test X_test, save y_test y_test 

    %     % non var, so axtime with hilbert on the transformed because it is
    %     % actually already filtered before transformation
    %     % matlab's hilbert
    %     hilbert_dat=[];
    %     for h=1:size(datatrn.trial,1) %trls
    %         for j=1:size(datatrn.trial,2) %ch
    %             hilbert_dat(h,j,:) = hilbert( squeeze(datatrn.trial(h,j,:)) ); % time should be in the first dimension.
    %         end
    %     end % power: is the squared magnitude of the complex vector at an instant in time. ... abs(hilbert_dat).^2 ... phase: angle(hilbert_dat)
    %     datatrn.trial = abs(hilbert_dat).^2;
    %     hilbert_dat=[];
    %     for h=1:size(datatst.trial,1) %trls
    %         for j=1:size(datatst.trial,2) %ch
    %             hilbert_dat(h,j,:) = hilbert( squeeze(datatst.trial(h,j,:)) ); % time should be in the first dimension.
    %         end
    %     end % power: is the squared magnitude of the complex vector at an instant in time. ... abs(hilbert_dat).^2 ... phase: angle(hilbert_dat)
    %     datatst.trial = abs(hilbert_dat).^2;
    %     hilbert_dat=[];
    %     for h=1:size(data.baseline,1) %trls
    %         for j=1:size(data.baseline,2) %ch
    %             hilbert_dat(h,j,:) = hilbert( squeeze(data.baseline(h,j,:)) ); % time should be in the first dimension.
    %         end
    %     end % power: is the squared magnitude of the complex vector at an instant in time. ... abs(hilbert_dat).^2 ... phase: angle(hilbert_dat)
    %     data.baseline = abs(hilbert_dat).^2;
    %     % performing baseline correction
    %     baseline_mat(:,:,1) = mean(data.baseline,3);
    %     baseline_ready = repmat(baseline_mat,1,1,size(data.trial,3));
    %
    %     datatrn.trial = ((datatrn.trial - baseline_ready(idx_trn,:,:))./baseline_ready(idx_trn,:,:))*100; % percentage change
    %     datatst.trial = ((datatst.trial - baseline_ready(idx_tst,:,:))./baseline_ready(idx_tst,:,:))*100;
    %
    %     % classification axtime time
    %     datatrn.trial = zscore(datatrn.trial, [], 1); datatst.trial = zscore(datatst.trial, [], 1);
    %     cfg=[];
    %     cfg.method = 'axtime'; cfg.classifier_type = {'lda'};
    %     cfg.do_parallel = 0;
    %     cfg.trn = datatrn; cfg.trn.trialinfo=datatrn.trialinfo';  cfg.folds=nan;
    %     cfg.tst = datatst; cfg.tst.trialinfo=datatst.trialinfo';
    %     auc(nn,:) = lv_classify(cfg);
    %

    %% scattering the data
    if nn==4
        % after CSP
        figure,
        scatter([datatrn.trial(datatrn.trialinfo==1,1);datatst.trial(datatst.trialinfo==1,1)], [datatrn.trial(datatrn.trialinfo==1,2);datatst.trial(datatst.trialinfo==1,2)]),
        hold on, scatter([datatrn.trial(datatrn.trialinfo==2,1);datatst.trial(datatst.trialinfo==2,1)], [datatrn.trial(datatrn.trialinfo==2,2);datatst.trial(datatst.trialinfo==2,2)],'r'),
        %         % var without CSP
        %                 figure,
        %                 scatter([hold_trn(datatrn.trialinfo==1,1);hold_tst(datatst.trialinfo==1,1)], [hold_trn(datatrn.trialinfo==1,2);hold_tst(datatst.trialinfo==1,2)]),
        %                 hold on, scatter([hold_trn(datatrn.trialinfo==2,1);hold_tst(datatst.trialinfo==2,1)], [hold_trn(datatrn.trialinfo==2,2);hold_tst(datatst.trialinfo==2,2)],'r'),
        %                 legend({'right hand', 'left hand'})
        %                 xlabel('C3'), ylabel('C4')
    end
auc
end

% group lvl TF
TF_struct.trial = TF_temp;
TF_struct.parametric = 1;
lv_tf(TF_struct, 1, 1); %data, do_stats, do_plot


x = randn(10,10);
[V,D] = eig(x)
[U,S,V] = svd(x);



lv_pretty_errorbar(data.time,auc,(auc.*0)+0.5, 2);


lv_pretty_errorbar( auc,(auc.*0)+0.5 );


lv_pretty_errorbar(data.time,left_hemi_c1,left_hemi_c2, 2); figure,
lv_pretty_errorbar(data.time,right_hemi_c1,right_hemi_c2, 2);



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