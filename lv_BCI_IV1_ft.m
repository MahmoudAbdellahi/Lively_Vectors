% lv_BCI_IV1_ft coverts BCI competition IV1 dataset to fieldtrip format
% that is ready for cleaning and other analyses
clear all
%load data and format
sbj_names = {'b','c','d','e','g'}; % removed a and f because not left vs right hand
raw_path = 'D:\sul''s code\Matt\sleep\erps\Organised\New exp\BCI_IV\Raw and filtered\BCICIV_calib_ds1';


%%
for nn=1:length(sbj_names)
    load([raw_path sbj_names{nn} '.mat']);
    dat.label = nfo.clab';
    dat.trial{1,1} = 0.1*double(cnt');
    dat.sampleinfo = [1 size(dat.trial{1,1},2)];
    dat.fsample = nfo.fs;
    dat.time{1,1} = 0:1/dat.fsample:(dat.sampleinfo(2)-1)/dat.fsample;
    
    %% filtering the continuous data
    cfg=[]; %the 100hz version was already band pass filtered and band stop
    %filtered
    %cfg.bpfilter  = 'yes'; because the data was already filtered 0.05 and 200 Hz
    %cfg.bpfreq    = [0.1 50]; % for new data sets consider FIR because it's better for offline analyses and not data driven as butter and IIR !
    %cfg.bsfilter  = 'yes'; % band-stop method
    %cfg.bsfreq    = [48 52];
    cfg.demean      = 'yes';
    dat  = ft_preprocessing(cfg, dat);
    
    %% segmenting into trials
    % because they rested sometimes between trials we will take the onset
    % of the cue and segment based on that ... because we extend 8sec.
    % after an onset with only the first 4sec. containing the activity
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
    data.trialinfo(data.trialinfo==-1)=2; % to make classes 1 and 2 was (-1 for class one or 1 for class two) so left hand is actually 2 now    % saving for cleaning and other analyses
    %clear cnt mrk nfo dat;
    %     lv_save([raw_path sbj_names{nn} '_filtered'], data, 'trial');
    % baseline
    baseline=3; onsets = mrk.pos'; onsets = onsets-baseline*dat.fsample; % for baseline we go back a bit
    onsets = repmat(onsets,1,baseline*dat.fsample+1);
    onsets = onsets + repmat(0:(baseline*dat.fsample),size(onsets,1),1);
    for i=1:size(dat.trial{1,1},1) % looping on channels and putting the trlxtime matrix
        temp= dat.trial{1,1}(i,:);
        data.baseline(:,i,:) = temp(onsets);
    end
    %     %% ERP analysis just to see if the effect is really lateralised
    % %         left_hemi_c1 = median(data.trial(data.trialinfo==2,ismember(data.label,'C3'),:), 1);
    % %         right_hemi_c1 = median(data.trial(data.trialinfo==2,ismember(data.label,'C4'),:), 1);
    % %         left_hemi_c2 = median(data.trial(data.trialinfo==1,ismember(data.label,'C3'),:), 1);
    % %         right_hemi_c2 = median(data.trial(data.trialinfo==1,ismember(data.label,'C4'),:), 1);
    % %         figure,
    % %         subplot(121), plot(data.time, [squeeze(left_hemi_c1)';squeeze(left_hemi_c2)']);
    % %         subplot(122), plot(data.time, [squeeze(right_hemi_c1)';squeeze(right_hemi_c2)']);
    % %         legend(nfo.classes) % not lateralised the effect is really in the variance/power
    % %         continue;
    %% scattering variances
    % %         cfg_preprocessing = [];  % we should do that on the continuous but just to see the effect
    % %         cfg_preprocessing.bpfilter = 'yes';
    % %         cfg_preprocessing.bpfilttype = 'fir';
    % %         cfg_preprocessing.bpfreq = [8 13];
    % %         cycles = 3; order = (cycles/min(cfg_preprocessing.bpfreq))*dat.fsample;
    % %         cfg_preprocessing.bpfiltord = round(order);
    % %         result= ft_preprocessing(cfg_preprocessing, data);
    % %
    % %         left_hemi_c1 = var(result.trial(data.trialinfo==2,ismember(data.label,'C3'),:),[], 3);
    % %         right_hemi_c1 = var(result.trial(data.trialinfo==2,ismember(data.label,'C4'),:),[], 3);
    % %         left_hemi_c2 = var(result.trial(data.trialinfo==1,ismember(data.label,'C3'),:),[], 3);
    % %         right_hemi_c2 = var(result.trial(data.trialinfo==1,ismember(data.label,'C4'),:),[], 3);
    % %         figure,
    % %         scatter(log(left_hemi_c1),log(right_hemi_c1))%high_low because ERD
    % %         hold on,
    % %         scatter(log(left_hemi_c2),log(right_hemi_c2))%low_high
    % %         legend(nfo.classes), xlabel('left channel'), ylabel('right channel')
    % %         continue;
    %% Hilbert transform
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
    figure,
    right_hemi_c1(nn,:) = median(data.trial(data.trialinfo==2,ismember(data.label,'C4'),:), 1); % to see power change from baseline axtime
    left_hemi_c1(nn,:) = median(data.trial(data.trialinfo==2,ismember(data.label,'C3'),:), 1);
    
    right_hemi_c2(nn,:) = median(data.trial(data.trialinfo==1,ismember(data.label,'C4'),:), 1);
    left_hemi_c2(nn,:) = median(data.trial(data.trialinfo==1,ismember(data.label,'C3'),:), 1);
    figure,
    subplot(121), plot(data.time, [squeeze(left_hemi_c1);squeeze(left_hemi_c2)]);
    subplot(122), plot(data.time, [squeeze(right_hemi_c1);squeeze(right_hemi_c2)]);
    legend(nfo.classes) % not lateralised the effect is really in the variance/power
    
    % classification across time
    %         cfg=[];
    %         cfg.method = 'axtime'; cfg.classifier_type = {'lda'};
    %         cfg.do_parallel = 0;
    %         cfg.trn = data; cfg.trn.trialinfo=data.trialinfo';  cfg.folds=5;
    %         auc(nn,:) = lv_classify(cfg);
    continue;
    %% CSP and Classification
    % CSP is supervised so we need to prevent info. leakage and need to use
    % only the training set to build el csp filters so we will split the
    % data into two data sets training and testing ..
    [datatrn,datatst] = deal(data);
    rng(1),
    idx = randperm(size(data.trial,1));
    idx_trn = idx(1:round(size(data.trial,1)/2)); idx_tst = idx(round(size(data.trial,1)/2)+1:size(data.trial,1));
    
    datatrn.trial = data.trial(idx_trn,:,:); datatst.trial = data.trial(idx_tst,:,:);
    datatrn.trialinfo = data.trialinfo(idx_trn); datatst.trialinfo = data.trialinfo(idx_tst);
    
    % csp
    warning('be careful of edges'); % filtering
    cfg_preprocessing                 = [];
    cfg_preprocessing.bpfilter        = 'yes';
    cfg_preprocessing.bpfreq          = [8 13];
    data_bp= ft_preprocessing(cfg_preprocessing, datatrn);
    cfg=[]; % building CSP filters
    cfg.method='csp';
    cfg.step='calculate';
    cfg.data=data_bp;
    cfg.data.trialinfo=datatrn.trialinfo';
    comp = lv_component_analysis(cfg);
    comp.id = comp.id(1:2); % choosing csp filters from sorted
    
    cfg.step='transform'; % applying the chosen filters to transform data,, in a CV paradigm when we use all data
    % to build the csp filters there will be leakage but we don't care much about it here
    cfg.comp=comp;
    datatrn.trial = lv_component_analysis(cfg);
    
    % apply to test
    cfg_preprocessing                 = [];
    cfg_preprocessing.bpfilter        = 'yes';
    cfg_preprocessing.bpfreq          = [8 13];
    data_bp= ft_preprocessing(cfg_preprocessing, datatst);
    cfg.data=data_bp;
    cfg.step='transform'; % applying the chosen filters to transform data,, in a CV paradigm when we use all data
    % to build the csp filters there will be leakage but we don't care much about it here
    cfg.comp=comp;
    datatst.trial = lv_component_analysis(cfg);
    
    % apply to baseline
    baselinedat = data;  baselinedat.trial = data.baseline; baselinedat.baseline=[];
    baselinedat.time = -baseline:1/dat.fsample:0;
    cfg_preprocessing                 = [];
    cfg_preprocessing.bpfilter        = 'yes';
    cfg_preprocessing.bpfreq          = [8 13];
    data_bp= ft_preprocessing(cfg_preprocessing, baselinedat);
    cfg.data=data_bp;
    cfg.step='transform'; % applying the chosen filters to transform data,, in a CV paradigm when we use all data
    % to build the csp filters there will be leakage but we don't care much about it here
    cfg.comp=comp;
    data.baseline = lv_component_analysis(cfg);
    
    
    
    % var or non var
    % var
    % %     datatrn.trial = log(var(datatrn.trial,[],3));
    % %     datatst.trial = log(var(datatst.trial,[],3));
    % %     % classification at a time
    % %     datatrn.trial = zscore(datatrn.trial, [], 1); datatst.trial = zscore(datatst.trial, [], 1);
    % %     cfg=[];
    % %     cfg.method = 'axtime'; cfg.classifier_type = {'lda'};
    % %     cfg.do_parallel = 0;
    % %     cfg.trn = datatrn; cfg.trn.trialinfo=datatrn.trialinfo';  cfg.folds=nan;
    % %     cfg.tst = datatst; cfg.tst.trialinfo=datatst.trialinfo';
    % %     auc(nn,:) = lv_classify(cfg);
    % non var, so axtime with hilbert on the transformed because it is
    % actually already filtered before transformation
    % matlab's hilbert
    hilbert_dat=[];
    for h=1:size(datatrn.trial,1) %trls
        for j=1:size(datatrn.trial,2) %ch
            hilbert_dat(h,j,:) = hilbert( squeeze(datatrn.trial(h,j,:)) ); % time should be in the first dimension.
        end
    end % power: is the squared magnitude of the complex vector at an instant in time. ... abs(hilbert_dat).^2 ... phase: angle(hilbert_dat)
    datatrn.trial = abs(hilbert_dat).^2;
    hilbert_dat=[];
    for h=1:size(datatst.trial,1) %trls
        for j=1:size(datatst.trial,2) %ch
            hilbert_dat(h,j,:) = hilbert( squeeze(datatst.trial(h,j,:)) ); % time should be in the first dimension.
        end
    end % power: is the squared magnitude of the complex vector at an instant in time. ... abs(hilbert_dat).^2 ... phase: angle(hilbert_dat)
    datatst.trial = abs(hilbert_dat).^2;
    hilbert_dat=[];
    for h=1:size(data.baseline,1) %trls
        for j=1:size(data.baseline,2) %ch
            hilbert_dat(h,j,:) = hilbert( squeeze(data.baseline(h,j,:)) ); % time should be in the first dimension.
        end
    end % power: is the squared magnitude of the complex vector at an instant in time. ... abs(hilbert_dat).^2 ... phase: angle(hilbert_dat)
    data.baseline = abs(hilbert_dat).^2;
    % performing baseline correction
    baseline_mat(:,:,1) = mean(data.baseline,3);
    baseline_ready = repmat(baseline_mat,1,1,size(data.trial,3));
    
    datatrn.trial = ((datatrn.trial - baseline_ready(idx_trn,:,:))./baseline_ready(idx_trn,:,:))*100; % percentage change
    datatst.trial = ((datatst.trial - baseline_ready(idx_tst,:,:))./baseline_ready(idx_tst,:,:))*100;
    
    % classification axtime time
    datatrn.trial = zscore(datatrn.trial, [], 1); datatst.trial = zscore(datatst.trial, [], 1);
    cfg=[];
    cfg.method = 'axtime'; cfg.classifier_type = {'lda'};
    cfg.do_parallel = 0;
    cfg.trn = datatrn; cfg.trn.trialinfo=datatrn.trialinfo';  cfg.folds=nan;
    cfg.tst = datatst; cfg.tst.trialinfo=datatst.trialinfo';
    auc(nn,:) = lv_classify(cfg);
    
    %% RNN
    
    %% Riemannian
    
    %% CORAL
    
    %% Euclidean alignment
    
end

lv_pretty_errorbar(data.time,auc,(auc.*0)+0.5, 2);


lv_pretty_errorbar( auc,(auc.*0)+0.5 );


lv_pretty_errorbar(data.time,left_hemi_c1,left_hemi_c2, 2); figure,
lv_pretty_errorbar(data.time,right_hemi_c1,right_hemi_c2, 2);

