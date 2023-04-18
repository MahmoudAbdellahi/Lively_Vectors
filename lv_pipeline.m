%% main .. calls other functions
% general params
start_up
addDirs

type = 'sleep'; % 'sleep' 'img'
% naming convention: dataset_convention = 'part*_sleep?'; hypno_convention = 'psgHypno-part*_sleep?'; * for participant no. and ? for part no.
sleep_stage = 3;
if strcmp(type,'img')==1, rawdir = 'E:\MRI_jittered_exp\allIMG'; sleep_stage=0; end % imagery doesn't have parts if your data doesn't have parts call it imagery
stats_window = [-4 4]; % this is the time window for which we do the stats. for determining outliers and interpolation

cleaned_path = ['D:\sul''s code\Matt\sleep\erps\Organised\New exp\RawData\allSleep\data\cleaned\final_cleaned_after_inspection\' mri_append 'part'];
preprocdir = 'D:\sul''s code\Matt\sleep\erps\Organised\New exp\RawData\allSleep\data\cleaned\';

fstruct = dir([rawdir '/' mri_append 'part*.eeg']); h = struct2cell(fstruct);
sbj = unique(str2double( (cellfun(@(x) (strtok(x,[mri_append 'part_'])), (h(1,:))','Un', 0))' ));

if isempty(mri_append), sbj(ismember(sbj,[ 2 19 13 25]))=[]; else sbj(ismember(sbj,[10 25]))=[]; end% bad sbj to remove from before MRI data,
% sbj2 had different reference and channels.. 19 13 25


% those are the ppnts that got sleep and imagery data in MRI data
sbj_MRI = [2 3 5 6 7 8 9 11 13 14 15 16 17 18 19 20 21 22 23 24  26 27 28 29 30 31 32 33];

% aggregating before MRI and MRI ppnts
sbj = [sbj sbj_MRI];
%% segmenting and cleaning
for nn=1:numel(sbj)
    data_parts = lv_check_parts(sbj(nn),type);
    cleaned_data = cell(1,data_parts);
    for part=1:data_parts
        cleaned_data{1,part} = lv_segment_filter_raw(sbj(nn),type, part, sleep_stage);
    end

    cleaned_data = cleaned_data(cell2mat(cellfun(@(x) (~isempty(x)),cleaned_data,'Un',0)));
    cleaned_data = ft_appenddata([],  cleaned_data{:}); % appending all parts together

    cleaned_data = lv_clean_segmented(cleaned_data, [], stats_window, sbj(nn));

    % convert to h5 for faster execution
    %save([preprocdir 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], 'cleaned_data', '-v7.3');

    lv_save([preprocdir mri_append 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], cleaned_data, 'trial');
end


%% manual artifact rejection for trials and channels
% if we have different data and all is h5 we will need to use lv_save/load
% instead of the save/load
for nn=1:numel(sbj)
    fprintf(['\n Manual artifact rejection for trials and channels, subject: ' num2str(sbj(nn)) '\n']);
    data = lv_load([preprocdir mri_append 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], 'trial');
    %load([preprocdir 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], 'cleaned_data');
    %data = cleaned_data; clear cleaned_data; % data is a variable passed to the lv_manual_cleaning script not function to debug the code manually because it is manual inspection

    % fixing IDs (if needed).. to make them unique because they repeat because of
    % having different parts .. will encode it as real the trls and imaginary the part
    if size(data.trialinfo,2)==4
        part_shift = [0 find(diff(data.trialinfo(:,2))<0) length(data.trialinfo(:,2))]; no_prts=1:length(part_shift)-1;
        if lv_check_parts(sbj(nn),type)~=length(no_prts), error('lv: parts mismatch between what is in data and what is on hard..'); end
        for j=1:length(no_prts), temp_idx=part_shift(j)+1:part_shift(j+1); % to start from idx 1
            data.trialinfo(temp_idx,5)= j; end % col. 5 for part no.
    end

    lv_manual_cleaning
    save([preprocdir 'possible_bad_channels\' mri_append 'part' num2str(sbj(nn)) '_' type '_possible_bad_ch_N' num2str(sleep_stage)], 'possible_bad_ch', '-v7.3');
    save([preprocdir 'manual_inspection_bad_trials\' mri_append 'part' num2str(sbj(nn)) '_' type '_manual_bad_trls_N' num2str(sleep_stage)], 'unique_id_bad_trls', '-v7.3');

    remaining_trls(nn,1)=size(data.trial,1)-size(unique_id_bad_trls,1);
    fprintf(['\n remaining_trls, subject: ' num2str(sbj(nn)) ', ' num2str(remaining_trls(nn,1)) '\n']);
    %save([preprocdir 'part' num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)], 'data', '-v7.3'); % cleaned data for further analyses
    fprintf(['\n Done. \n']);
end


%% loading the saved and mixing with martyna's cleaning
for nn=1:numel(sbj)
    fprintf(['\n Manual artifact rejection MIXING, subject: ' num2str(sbj(nn)) '\n']);
    %     load([preprocdir mri_append 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], 'cleaned_data'); % old save/load left commented because sleep data was saved this way
    %     data = cleaned_data; clear cleaned_data; % data is a variable passed to the lv_manual_cleaning script not function to debug the code manually because it is manual inspection
    data = lv_load([preprocdir mri_append 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], 'trial');

    % fixing IDs .. to make them unique because they repeat because of
    % having different parts .. will encode it as real the trls and imaginary the part
    if size(data.trialinfo,2)==4
        part_shift = [0 find(diff(data.trialinfo(:,2))<0) length(data.trialinfo(:,2))]; no_prts=1:length(part_shift)-1;
        if lv_check_parts(sbj(nn),type)~=length(no_prts), error('lv: parts mismatch between what is in data and what is on hard..'); end
        for j=1:length(no_prts), temp_idx=part_shift(j)+1:part_shift(j+1); % to start from idx 1
            data.trialinfo(temp_idx,5)= j; end % col. 5 for part no.
    end
    if sleep_stage==0 % wake
        load([preprocdir 'manual_inspection_bad_trials\' mri_append 'part' num2str(sbj(nn)) '_' type '_manual_bad_trls_N' num2str(sleep_stage)], 'unique_id_bad_trls');
        mix_bad = unique_id_bad_trls;
    else
        % loading martyna's bad trials and saved bad trials
        %load([preprocdir 'manual_inspection_bad_trials\Martyna\data_cleaned\part' num2str(sbj(nn)) '_' type '_manual_bad_trls2_N' num2str(sleep_stage)], 'unique_id_bad_trls2');
        load([preprocdir 'manual_inspection_bad_trials\' mri_append 'part' num2str(sbj(nn)) '_' type '_manual_bad_trls_N' num2str(sleep_stage)], 'unique_id_bad_trls');
        mix_bad = unique_id_bad_trls; %union(unique_id_bad_trls2,unique_id_bad_trls,'rows'); % intersect/union applies unique automatically
    end
    remaining_trls(nn,1)=size(data.trial,1)-size(mix_bad,1);
    [~,bad_trls] = intersect([data.trialinfo(:,2) data.trialinfo(:,5)] , mix_bad, 'rows'); % idx in trialinfo


    idx = 1:size(data.trial,1); idx(bad_trls)=[]; cfg = []; cfg.trials = idx;
    data = ft_selectdata(cfg, data);

    fprintf(['\n remaining_trls, subject: ' num2str(sbj(nn)) ', ' num2str(remaining_trls(nn,1)) '\n']);

    lv_save([preprocdir 'final_cleaned_after_inspection\' mri_append 'part' num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)], data, 'trial'); % cleaned data for further analyses
    fprintf(['\n Done. \n']);
end
remaining_trls


%% if cogent was used make sure to use this block
% remove possible hypothetical triggers sent by cogent in the beginning of recodring of parts
% possible removal of the first trial of each part because cogent was sending a hypothetical trigger in the beginning
% check the record by eye and then change this variable to the trial(s) to reject and move in this code line by line
% % for nn=9:numel(sbj)
% %     sbj(nn)
% %     cleaned_data = lv_load([preprocdir 'final_cleaned_after_inspection\' mri_append 'part' num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)],'trial');
% %     % inspect manually to see if the first trial of the part (record) is labeled 1 ..
% %     part=unique(cleaned_data.trialinfo(:,end)); record=[];
% %     for i=1:length(part), record(i)= min(find(cleaned_data.trialinfo(:,end)== part(i))); end
% %     remove = [155]; % check the record by eye and then change this variable to the trial(s) to reject and move in this code line by line
% %     idx=1:size(cleaned_data.trial,1); idx(remove)=[];cfg = []; cfg.trials = idx;
% %     data = ft_selectdata(cfg, cleaned_data);
% %
% %     sbj(nn)
% %     % if you have trials to remove and the data is updated then save
% %     lv_save([preprocdir 'final_cleaned_after_inspection\' mri_append 'part' num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)], data, 'trial'); % cleaned data for further analyses
% %     fprintf(['\n Done. \n']);
% % end

%% also do the exact segmenting in lv_segment raw and remove any trial with short duration that doesn't make sense
% because the triggers are sometimes incorrectly sent very abruptly

%% ERP analysis
erps_temp=[];
for nn=1:numel(sbj)
    fprintf(['\n ERP analysis, subject: ' num2str(sbj(nn)) '\n']);
    cleaned_data = lv_load([cleaned_path num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)],'trial');

    cfg=[]; cfg.latency=[-0.5 1.5]; cleaned_data=ft_selectdata(cfg,cleaned_data);% reducing to time of trial
    % putting trials of left hand together and for right hand as well
    classes = unique(cleaned_data.trialinfo(:,1)); cleaned_data.trialinfo(ismember(cleaned_data.trialinfo(:,1),classes(1:2)),1 ) = 1; cleaned_data.trialinfo(ismember(cleaned_data.trialinfo(:,1),classes(3:4)),1 ) = 2;

    erps = lv_erp(cleaned_data, 0, 0); %data, do_stats, do_plot ... returns 2_ch_time the first rpt is cond1 erp then cond2
    erps_temp = [erps_temp ; erps.trial]; % erps_temp aggregates all the erps of different sbj
end

% group lvl ERP
erps.trial = erps_temp;
lv_erp(erps, 0, 1); %data, do_stats, do_plot


%% TF analysis
TF_temp=[];
for nn=1:numel(sbj)
    fprintf(['\n TF analysis, subject: ' num2str(sbj(nn)) '\n']);
    cleaned_data = lv_load([cleaned_path num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)],'trial');

    cfg=[]; cfg.latency=[-0.5 1.5]; cleaned_data=ft_selectdata(cfg,cleaned_data);% reducing to time of trial
    % putting trials of left hand together and for right hand as well
    classes = unique(cleaned_data.trialinfo(:,1)); cleaned_data.trialinfo(ismember(cleaned_data.trialinfo(:,1),classes(1:2)),1 ) = 1; cleaned_data.trialinfo(ismember(cleaned_data.trialinfo(:,1),classes(3:4)),1 ) = 2;

    [ TF_struct ] = lv_tf(cleaned_data, 0, 0); %data, do_stats, do_plot .. gets the TF in TF_struct.trial
    TF_temp = [TF_temp ; TF_struct.trial]; % erps_temp aggregates all the erps of different sbj
end

TF_temp_bl = TF_temp;  % with baseline
save TF_temp_bl TF_temp_bl;
save TF_temp TF_temp % 22sbj TF analyses

% group lvl TF
TF_struct.trial = TF_temp;
lv_tf(TF_struct, 0, 1); %data, do_stats, do_plot


%% ITPC

% the reason is to include ERP, TF, and ITPC and be able to analyse the
% features to see where differences between classes exist without
% classification

%% feature extraction and classification
%% running a classification pipeline
% sbj = [3 4 5]; % just to test .. when you actually do this use the main cell in this file the first one to load the params

trn_ds = 'img';
tst_ds = 'sleep';
if strcmp(trn_ds,'img')==1, sleep_stage_trn=0; else trn_ds='sleep'; sleep_stage_trn=3; end
if strcmp(tst_ds,'img')==1, sleep_stage_tst=0; else tst_ds='sleep'; sleep_stage_tst=3; end
sequence = {'lv_feature_extractor','erp 100'}; % to run many sequences back to back put new seq. in rows so it will be run(s) x sequence
[trn,tst,classification_result,fidelity_pt,fidelity_pt2,lens,lens_sleep,recurrence_distrib]=deal([]);
temp=[]; erp=[]; erp_trn=[];  spindle_likelihood=[];
vschance = 1; % vs chance or do we have control condition
if vschance==1, conditions=1; else, conditions=2; end
cleaned_path = ['D:\sul''s code\Matt\sleep\erps\Organised\New exp\RawData\allSleep\data\cleaned\final_cleaned_after_inspection\part'];
if ~exist('sequence','var') % if the sequence not provided then load from excel sheet
    table_specs = readtable([pwd '\lv_automation_runs.xlsx']);
    titles = table_specs.Properties.VariableNames;
    variables = table_specs.Variables;
    run_id=[];
    for i=1:size(variables,1), if strcmp(char(variables(i,end)),'done')==0, run_id=[run_id ; i]; end, end
    sequence = variables(run_id,2:end);
    fprintf(['\n lv: Applying pipeline sequence(s): \n']); disp(sequence);
    record=repmat({' '},size(variables,1),1); record(run_id,1)={'done'};
    table_specs = [table_specs record];
    writetable(table_specs, [pwd '\lv_automation_runs.xlsx']); % will put done here better than checking again at the end of execution
end
load lv_layout lv_layout; axtimeAcc=[];
pathname = 'D:\sul''s code\Matt\sleep\erps\Organised\New exp'; % path to save the figures of results
for run=1:size(sequence,1)
    for cond=1:conditions % example: exp, adp
        for nn=1:numel(sbj)
            if nn>=22 % because the before MRI data got 21 ppnts
                % load MRI data and continue looping
                cleaned_path = ['D:\sul''s code\Matt\sleep\erps\Organised\New exp\RawData\allSleep\data\cleaned\' ...
                    'final_cleaned_after_inspection\MRI_part'];
            end
            fprintf(['\n lv: Working on pipeline for sbj: ' num2str(sbj(nn)) ' \n']);
            trn.data = lv_load([cleaned_path num2str(sbj(nn)) '_' trn_ds '_manual_cleaned_N' num2str(sleep_stage_trn)],'trial');
            tst.data = lv_load([cleaned_path num2str(sbj(nn)) '_' tst_ds '_manual_cleaned_N' num2str(sleep_stage_tst)],'trial');

            % cutting the trials with respect to their variable lengths
            temp = tst.data.trial_lens;
            [~,~,common_all] = intersect([tst.data.trialinfo(:,2) tst.data.trialinfo(:,5)],[temp(:,2) temp(:,5)], 'rows','stable');
            temp = temp(common_all,:);
            if ~isequal(tst.data.trialinfo(:,1:5),temp(:,1:5)), error('lv_ mismatch between records'); end
            tst.data.trialinfo = [tst.data.trialinfo temp(:,6)];
            % rejecting bad trials and fixing long ones
            temp = tst.data.trialinfo(:,6);
            bad_id = find(temp<2.5);
            if length(bad_id)>0, bad_id, warning('trial <2.5'); end % not real error, but just to see if this happens
            tst.data.trial(bad_id,:,:)=[]; tst.data.trialinfo(bad_id,:)=[]; tst.data.sampleinfo(bad_id,:)=[];
            temp = tst.data.trialinfo(:,6);
            long_id = find(temp>3.5);
            tst.data.trialinfo(long_id,6)=3.5; clear temp;

            % %             % baseline correction
            % %             tst_bl.data = tst.data; bl=[-4 0];
            % %             tst_bl.data.trial = bsxfun(@minus,tst_bl.data.trial, mean(tst_bl.data.trial(:,:,...
            % %                 nearest(tst.data.time,bl(1)):nearest(tst.data.time,bl(2))), 3) );
            % %             tst.data = tst_bl.data;
            % %
            % %             bl=[-4 0];
            % %             trn.data.trial = bsxfun(@minus,trn.data.trial, mean(trn.data.trial(:,:,...
            % %                 nearest(trn.data.time,bl(1)):nearest(trn.data.time,bl(2))), 3) );

            %             % Spindle and SO prediction
            %             cfg=[]; cfg.latency=[-3.04 3.04]; dat{1,1} = ft_selectdata(cfg, tst.data); spindle=1;
            %             result = lv_SO_Spindle_FeatureExtraction(dat, spindle); mx_dataClick_1 = result;
            %             if spindle==0, ch_id = find(ismember(tst.data.label,'fz'));
            %                 X = mx_dataClick_1{1,1}{ch_id, 10}.features;
            %                 Z = X( ~isnan( mean(X,2) ) , : ); trl_id = 1:size(tst.data.trial,1); trl_id = trl_id(~isnan( mean(X,2) ))';
            %                 mx_dataClick_1{ch_id, 10}.features = Z;
            %                 features_test{1,1} = [cos(mx_dataClick_1{ch_id, 10}.features(:,1))  sin(mx_dataClick_1{ch_id, 10}.features(:,1))   mx_dataClick_1{ch_id, 10}.features(:,2:end)];
            %             else
            %                 X = result'; features_test{1,1} = cell2mat(X(1,:)); trl_id = (1:size(tst.data.trial,1))';
            %             end
            %
            %             ids = lv_SO_Spindle_classification(features_test, trl_id, spindle);
            %             tst.data.trial = tst.data.trial(ids,:,:); tst.data.trialinfo = tst.data.trialinfo(ids,:); tst.data.sampleinfo = tst.data.sampleinfo(ids,:);
            %             no_trls(nn,1) = size(tst.data.trial,1);


            % time frequency analysis
%             cfg=[]; cfg.latency = [-0.6 3.8]; tst.data = ft_selectdata(cfg, tst.data);
%             TF_struct= do_tf(tst.data, [-0.3 -0.1], [1 30]);
% 
%             TF_struct.trial = squeeze(mean(TF_struct.trial, 1)); % average of channels
%             tf_pw(nn,:,:) = TF_struct.trial; % participant x frequency x time
% 
%             continue;

            % ERP analysis
            baseline_id = nearest(tst.data.time,-0.2):nearest(tst.data.time,0); % time 

            baseline = (mean(tst.data.trial(:,:,baseline_id),3));
            baseline = repmat(baseline, 1,1,size(tst.data.trial,3));
            erp(nn,:) = mean(mean(tst.data.trial - baseline, 2),1);

            continue;

            % extracting features
            cfg=[];
            cfg.data=trn.data;
            cfg.sequence = sequence;
            result_trn = lv_build_pipeline(cfg); % for trn

            cfg.data=tst.data;
            result_tst = lv_build_pipeline(cfg); % for tst

            % making classes binary, will not do this if data was different
            classes=unique(result_trn.trialinfo(:,1)), if length(classes)~=4, error('lv: found classes are not 4 !'); end% the first two are aggregated together to be left hand and then the second two as right hand
            result_trn.trialinfo( ismember(result_trn.trialinfo(:,1),classes([1])) ,1)  = 1;
            result_trn.trialinfo( ismember(result_trn.trialinfo(:,1),classes([2])) ,1)  = 2;
            result_trn.trialinfo( ismember(result_trn.trialinfo(:,1),classes([3])) ,1)  = 3;
            result_trn.trialinfo( ismember(result_trn.trialinfo(:,1),classes([4])) ,1)  = 4;

            classes=unique(result_tst.trialinfo(:,1)), if length(classes)~=4, error('lv: found classes are not 4 !'); end% the first two are aggregated together to be left hand and then the second two as right hand
            result_tst.trialinfo( ismember(result_tst.trialinfo(:,1),classes([1])) ,1)  = 1;
            result_tst.trialinfo( ismember(result_tst.trialinfo(:,1),classes([2])) ,1)  = 2;
            result_tst.trialinfo( ismember(result_tst.trialinfo(:,1),classes([3])) ,1)  = 3;
            result_tst.trialinfo( ismember(result_tst.trialinfo(:,1),classes([4])) ,1)  = 4;
            % data trimming
            cfg=[];
            %             cfg.channel = {'c1','c3','c5','cp3','fc3'...
            %                 'c6','c4','c2','cp4','fc4'}; % diamond shaped
            %             cfg.channel = lower({'C6','C4','C2','C1','C3','C5', 'CP5','CP3', 'CP1', 'CP2', 'CP4', 'CP6'...
            %                 'FC5','FC3','FC1', 'FC2','FC4','FC6'}); % aoi
            %                         cfg.channel = lower({'FC3','C5', 'C3', 'C1', 'CP3', 'FC4', 'C6', 'C4', 'C2', 'CP4'}); % martyna's
            % % %             hold = result_tst.trialinfo(:,1); hold2 = result_trn.trialinfo(:,1);
            % % %             cfg.channel = lower({'fp1','af3','f1','f3','f5','f7','fc1','fc3','fc5','ft7','c1','c3','c5','t7',...
            % % %                 'cp1','cp3','cp5','tp7','p1','p3','p5','p7','po3','po7','o1'});
            % % %             result_tst1=ft_selectdata(cfg,result_tst);
            % % %             result_trn1=ft_selectdata(cfg,result_trn);
            % % %
            % % %             cfg.channel = lower({'fp2','af4','f2','f4','f6','f8','fc2','fc4','fc6','ft8','c2','c4','c6','t8',...
            % % %                 'cp2','cp4','cp6','tp8','p2','p4','p6','p8','po4','po8','o2'});

            %             cfg.channel = lower({'fp2','af4','f2','f4','f6','f8','fc2','fc4','fc6','ft8','c2','c4','c6','t8',...
            %                 'cp2','cp4','cp6','tp8','p2','p4','p6','p8','po4','po8','o2'});
            % % %             result_tst=ft_selectdata(cfg,result_tst);
            % % %             result_trn=ft_selectdata(cfg,result_trn);
            % % %
            % % %             result_tst.trial = result_tst1.trial-result_tst2.trial;
            % % %             result_trn.trial = result_trn1.trial-result_trn2.trial;
            cfg=[];

            cfg.latency=[0 3.5];
            %result_trn=ft_selectdata(cfg,result_trn);
            result_tst=ft_selectdata(cfg,result_tst);
            if strcmp(trn_ds,'img')==1, cfg.latency=[0 1.1]; result_trn=ft_selectdata(cfg,result_trn);  end
            if strcmp(tst_ds,'img')==1, cfg.latency=[0 1.1]; result_tst=ft_selectdata(cfg,result_tst); end
            sz = size(result_tst.trial);
            % variable length nans
            idx = cell2mat(arrayfun(@(x) (nearest(result_tst.time,x)), result_tst.trialinfo(:,6),'Un',0));
            for i=1:size(result_tst.trial,1), result_tst.trial(i,:,idx(i)+1:end) = nan; end

            %             ratio=0.3333; % half of stims
            %             cfg=[]; h=round(size(result_trn.trial,1).* ratio);
            %             cfg.trials=1: h;
            %             result_trn = ft_selectdata(cfg,result_trn);

            % % % %             % crossfreq_PAC
            % % % %             for ch=1:length(tst.data.label)
            % % % %                 cfg=[]; cfg.channel=tst.data.label(ch); temp=ft_selectdata(cfg,tst.data);
            % % % %                 cfg=[];
            % % % %                 cfg.method = 'crossfreq_PAC'; cfg.data = temp;
            % % % %                 cfg.freqLow = [0.5 2]; cfg.freqHigh = [11 16];
            % % % %                 cfg.toi_window = [0 2.5]; % we cut to this length later and filter with longer to prevent edge artifacts
            % % % %                 result = lv_feature_extractor(cfg);
            % % % %                 pac_info(nn,ch,:) = [result.coupling_strength result.phase_of_coupling result.pval];
            % % % %             end
            % % % %             continue;


            % % %             result_tst.trialinfo(:,1) = hold; result_trn.trialinfo(:,1) = hold2;
            % fixing channels
            % % %             result_tst.lv_layout = lv_layout; result_tst.fix_option='auto';
            % % %             result_tst = lv_fix_channels(result_tst);
            % % %             result_trn.lv_layout = lv_layout;
            % % %             result_trn.fix_option='auto'; % this lowers the rank of data so that be carefu lto use ICA after that because now the interpoalted channel will be dependent on others
            % % %             result_trn = lv_fix_channels(result_trn);
            % % % %             if any( ismember(union(result_trn.suggested_bad_channel,result_tst.suggested_bad_channel),'cz'))
            % % % %                 continue; end % because this channel is used for spindle detection
            % % % % % %           Don't variance of channels is handeled by zscoring  bad_channels = union(result_trn.suggested_bad_channel,result_tst.suggested_bad_channel);
            % % % % % %             if ~isempty(bad_channels)
            % % % % % %                 cfg=[]; cfg.channel = result_trn.label(~ismember(result_trn.label,bad_channels));
            % % % % % %                 result_trn = ft_selectdata(cfg,result_trn); result_tst = ft_selectdata(cfg,result_tst);
            % % % % % %             end

            hold_trialinfo = result_tst.trialinfo;

            % %             % splitting trials based on the SO phase
            % %             cfg=[];
            % %             cfg.data = tst.data; cfg.phenom_ch = 'fz';
            % %             cfg.method = 'up_going';
            % %             stim_up = lv_reduce_trials(cfg);
            % %             result_tst.trial = result_tst.trial(stim_up==1,:,:);
            % %             result_tst.trialinfo = result_tst.trialinfo(stim_up==1,:);


            %             % phenomena detection
            %             %cfg=[]; cfg.latency=[0 2.5]; temp_data=ft_selectdata(cfg,tst.data);
            %             temp_data = tst.data;
            %             cfg=[];
            %             cfg.data = temp_data; cfg.phenom_ch = 'fz';
            %             cfg.method = 'up_going'; %'up_going';   ASC
            %             %cfg.lock_to_trough=1;
            %             %             cfg.lock_to_spindle_trough=1;
            %             %             cfg.lock_to_spindle=1;
            %             [ good_trials ] = lv_reduce_trials(cfg); %trials_kept(nn,1)=sum(good_trials(:));
            %
            %             % %             ids = 1:size(tst.data.trial,1); ids(corr_majority{nn,1})=[];
            %             % %             corr_majority{nn,1}=ids';
            %             % %             phase(nn,1)=mean(good_trials(corr_majority{nn,1}));  continue; % calculating the stimulation phase of specific trials
            %
            %             % cutting good_trials to match the timing of trials
            %             id = nearest(tst.data.time,0):nearest(tst.data.time,3.5);
            %             good_trials = good_trials(:,id);
            %             %             good_trials( sum(good_trials,2)>0 ,:) = 1;% trials splitting by making the trials with phenom all ones

            % %             % phenomena detection on all channels
            % %             % gets the SO phase and then gets the spindles on the upgoing
            % %             temp_data = tst.data; good_trials=[];
            % %             for i=1:size(tst.data.label)
            % %                 cfg=[]; % SO_phase
            % %                 cfg.data = temp_data; cfg.phenom_ch = tst.data.label(i);
            % %                 cfg.method = 'up_going';
            % %                 temp_so = lv_reduce_trials(cfg);
            % %                 cfg.method = 'ASC'; % Spindles
            % %                 cfg.lock_to_spindle=1;
            % %                 temp_spindles = lv_reduce_trials(cfg);
            % % %                 good_trials(:,i,:) = temp_so & temp_spindles; % trial ch time of flags
            % %                 temp_so_up=temp_so; temp_so_up(isnan(temp_so_up))=0;
            % %                 temp_so_down=temp_so; temp_so_down(temp_so_down==1)=0; temp_so_down(isnan(temp_so_down))=1;
            % %                 countOnUp=temp_so_up & temp_spindles; countOnUp=sum(countOnUp(:));
            % %                 countOnD=temp_so_down & temp_spindles; countOnD=sum(countOnD(:));
            % %                 spindle_likeli(nn,i) = sum(countOnUp(:))./(countOnUp+countOnD);% likelihood of spindles on the upgoing phase
            % %             end
            % %             continue;
            % %             good_trials = double(squeeze(sum(good_trials,2))>0); % or
            % %             % cutting good_trials to match the timing of trials
            % %             id = nearest(tst.data.time,0):nearest(tst.data.time,3.5);
            % %             good_trials = good_trials(:,id);
            % %             % cutting to good pts.
            % %             events_kept(nn,1) = sum(good_trials(:));
            % %             if events_kept(nn,1)<20, continue; end
            % %             % putting nans for the non phenom parts to exclude them
            % %             good_trials(good_trials==0)=nan;
            % %             temp(:,1,:)=good_trials; temp = repmat(temp,1,size(result_tst.trial,2),1);
            % %             result_tst.trial = temp .* result_tst.trial;



            %             keep_id = find(sum(good_trials,2)>0);
            %             spindle_likelihood(nn,:) = mean(good_trials,1);
            % time locking
            %             for i=1:size(good_trials,1)
            %                 temp = find(good_trials(i,:));
            %                 if length(temp)>1
            %                     good_trials(i,:) = zeros(1,length(good_trials(i,:))); good_trials(i,temp(1))=1;
            %                 end
            %             end
            % %             result_tst = lv_event_timelocked(result_tst, [0 0], good_trials);
            % %             result_tst.trial = permute(result_tst.trial,[1 3 2]);% in case of one time pt [0 0] timelocking
            % %             events_kept(nn,1)= size(result_tst.trial,1)
            % %             if events_kept(nn,1)<20, continue; end
            % non-time locking .. trials based on event
            %             result_tst.trial = result_tst.trial(keep_id,:,:);
            %             result_tst.trialinfo = result_tst.trialinfo(keep_id,:);
            %events_kept(nn,1)= size(result_tst.trial,1)
            %if events_kept(nn,1)<20, continue; end




            % % %             % filtering in spindle band
            % % %             cfg_preprocessing                 = [];
            % % %             cfg_preprocessing.bpfilter        = 'yes';
            % % %             cfg_preprocessing.bpfreq          = [11 16];
            % % %             data_bp= ft_preprocessing(cfg_preprocessing, tst.data  );
            % % %             cfg=[]; cfg.latency=[-2 0]; cfg.channel = lower({'FC3','C5', 'C3', 'C1', 'CP3', 'FC4', 'C6', 'C4', 'C2', 'CP4'}); % martyna's
            % % %             data_bp=ft_selectdata(cfg,data_bp); %1.860 2.4150
            % % %
            % % %             %data_bp.trial = data_bp.trial(:,ismember(result_tst.label,'cz'),:);
            % % %
            % % %             hilbert_dat=[];
            % % %             % matlab's hilbert
            % % %             for h=1:size(data_bp.trial,1) %trls
            % % %                 for j=1:size(data_bp.trial,2) %ch
            % % %                     hilbert_dat(h,j,:) = hilbert( squeeze(data_bp.trial(h,j,:)) ); % time should be in the first dimension.
            % % %                 end
            % % %             end % power: is the squared magnitude of the complex vector at an instant in time. ... abs(hilbert_dat).^2
            % % %             hilbert_dat = abs(hilbert_dat).^2; % trl x ch x time, matlab byfham en abs hwa magnitude el complex vector length el complex vector ya3ne
            % % %             D1 = squeeze(median(hilbert_dat,3));
            % % %             D1 = mean(D1,2);
            % % %             %[val,id] = sort(D1,'descend');
            % % %             %keep_id = id(1:round(0.1*length(D1)) );
            % % %             keep_id = find(D1 > median(D1,1));
            % % %
            % % % % % % % % %             h=round(length(D1)*0.3333); % three parts spindle power
            % % % % % % % % %             pw(nn,1)=mean(D1(1:h)); pw(nn,2)=mean(D1(h+1:round(length(D1)*0.6666)));
            % % % % % % % % %             pw(nn,3)=mean(D1((round(length(D1)*0.6666))+1:length(D1) )); continue;
            % % % % % %
            % % % % % % %             corr_majority{nn,1} = keep_id; continue;
            % % % % % %
            % % % % % % % % % % %             pw = squeeze(hilbert_dat); % splitting based on the power of pts not trials
            % % % % % % % % % % %             masked_pw = pw .* ((squeeze(result_tst.trial(:,1,:)).*0)+1); % masked pw
            % % % % % % % % % % %             masked_pw = masked_pw(:); masked_pw( masked_pw<nanmedian(masked_pw) )=nan;
            % % % % % % % % % % %             mask = reshape(masked_pw,size(pw)); temp(:,1,:)=mask; temp=repmat(temp,1,size(result_tst.trial,2),1);
            % % % % % % % % % % %             result_tst.trial = result_tst.trial .* temp;
            % % % % % %
            % % % % % %
            % % % % % %             % %             idx = 1:size(result_tst.trial,1); idx(keep_id)=[];
            % % %             result_tst.trial = result_tst.trial(keep_id,:,:);
            % % %             result_tst.trialinfo = result_tst.trialinfo(keep_id,:);


            % %             % PAC of SO and spindles keep_id, check page 207 Notes 9 important
            % %             tst.data.trial = tst.data.trial(keep_id,:,:);
            % %             tst.data.trialinfo = tst.data.trialinfo(keep_id,:); hold_trialinfo = result_tst.trialinfo; sz = size(result_tst.trial);
            % %             idx = cell2mat(arrayfun(@(x) (nearest(tst.data.time,x)), result_tst.trialinfo(:,6),'Un',0)); % variable length nans
            % %
            % %             keepidx = single(ones(size(tst.data.trial))); time0=nearest(tst.data.time,0);
            % %             for i=1:size(tst.data.trial,1)
            % %                 keepidx(i,:,idx(i)+1:end) = nan; keepidx(i,:,1:time0-1) = nan;
            % %             end
            % %             cfg=[]; cfg.data = tst.data;
            % %             cfg.method = 'PAC'; cfg.keepidx = keepidx; % this is for the parts out of edges and nans of jittering
            % %             [ result ] = lv_reduce_trials(cfg);
            % %             phase(nn,1) = result.phase; uniformity(nn,1) = result.uniformity;
            % % %             uniformity(nn,1) = result.coupling_strength;
            % % %             phase(nn,1) = result.phase_of_coupling;
            % % %             continue;

            % %             events_kept(nn,1)= size(result_tst.trial,1)
            % %             if events_kept(nn,1)<20, continue; end

            %             % surface laplacian high pass spatial filtering
            %             if nn==1, cfg=[]; lv_layout.pos = [lv_layout.pos zeros(length(lv_layout.label),1)]; end % dropping z dimension
            %             cfg.data=result_trn; cfg.method = 'spline'; cfg.elec = lv_layout; cfg.elec.elecpos = cfg.elec.pos;
            %             result_trn = ft_scalpcurrentdensity(cfg, result_trn);
            %             result_tst = ft_scalpcurrentdensity(cfg, result_tst);

            % checking topos
            %             cfg=[]; % converting time pt to topo.
            %             cfg.data=result_trn;
            %             cfg.method = 'eeg_to_topos_video';
            %             cfg.data.layout = lv_layout;
            %             peak_topo = lv_feature_extractor(cfg);


            % % %             % pooled ICA
            % % %             % ICA
            % % %             temp = [lv_3dto2d(result_trn.trial,size(result_trn.trial)) lv_3dto2d(result_tst.trial,size(result_tst.trial))];
            % % %             result = ft_appenddata([], result_trn,result_tst);
            % % %             comp_cfg = []; comp_cfg.numcomponent=rank(temp)
            % % %             comp_cfg.method = 'runica'; rng(1); clear temp;
            % % %             comp_cfg.demean = 'no';
            % % %             comp = ft_componentanalysis(comp_cfg, result); % components are calculated using trn only ..
            % % %             cfg=[]; cfg.channel = 1:10; comp = ft_selectdata(cfg,comp);
            % % %
            % % %             % project back
            % % % %             cfg_comp=[]; cfg_comp.component=[5:comp_cfg.numcomponent]; cfg_comp.demean = 'no';
            % % % %             result_trn = ft_rejectcomponent(cfg_comp, comp, result_trn); % projection on data
            % % % %             result_tst = ft_rejectcomponent(cfg_comp, comp, result_tst); % projection on data
            % % %
            % % %             trn_lbl = result_trn.trialinfo; tst_lbl = result_tst.trialinfo;
            % % %             cfg=[]; cfg.trials = 1:size(result_trn.trial,1);
            % % %             result_trn = ft_selectdata(cfg,comp); result_trn.trialinfo=trn_lbl;
            % % %             result_trn.trial = lv_2dto3d(cell2mat(result_trn.trial),[length(result_trn.trial) size(result_trn.trial{1,1},1) size(result_trn.trial{1,1},2)]);
            % % %
            % % %             cfg=[]; cfg.trials = size(trn_lbl,1)+1:length(comp.trial);
            % % %             result_tst = ft_selectdata(cfg,comp); result_tst.trialinfo=tst_lbl;
            % % %             result_tst.trial = lv_2dto3d(cell2mat(result_tst.trial),[length(result_tst.trial) size(result_tst.trial{1,1},1) size(result_tst.trial{1,1},2)]);

            % all pts position encoding
            %             cfg=[]; cfg.data=result_trn; cfg.data.layout=lv_layout; cfg.method = 'encoded_pts';
            %             result_trn = lv_feature_extractor(cfg); % not working


            result_trn.trial = single(result_trn.trial); result_tst.trial = single(result_tst.trial);
            %             result_trn.trial = zscore(result_trn.trial,[],1); result_tst.trial = zscore(result_tst.trial,[],1); % zscoring
            % spatiotemporal PCA .. no centering because PCA data is
            % demeaned also PCA need the cov matrix and to get the cov data is centered so the data when we get
            % the PCs is centered and when we apply because it was demeaned
            % claculate
            cfg=[];
            cfg.data = result_tst;
            %cfg.data2 = result_trn;
            cfg.method = 'pca';
            cfg.step = 'calculate';
            cfg.centered = 1; % data is centered with demean in preprocessing
            comp = lv_component_analysis(cfg); %cfg=rmfield(cfg,'data2');
            % transform
            cfg.eigVects = comp.eigVects;
            cfg.chosen = zeros(length(comp.eigVals),1);
            %             cfg.chosen(1:4) = 1;
            id = find(cumsum(comp.eigVals)> 0.95); % 80% of variance
            cfg.chosen(1:id(1))=1;

            % %             eigVectors(:,nn) = cfg.eigVects(:,1);
            % %             continue;
            %             cfg.chosen(comp.eigVals>1)=1; % Kaiser and change
            %             the lv_component to return the actual eigVal not the explained variance and zscore before doing this
            cfg.step = 'transform';
            result_tst.trial = lv_component_analysis(cfg);
            cfg.data = result_trn;
            result_trn.trial = lv_component_analysis(cfg);


            %             % temporal compression with wilson's method
            %             cfg=[]; cfg.data_long=result_tst; cfg.data_short=result_trn; cfg.method='compression_wilson_jittered';
            %             cfg.centered = 0; % for PCA because it will need to center each ch_time
            %             warning('for compression analysis no smoothing should be done and the length of trials should be the same');
            %             temp = lv_post_classification(cfg);
            %             compressions{1,nn} = temp.compressionRatio(temp.score==max(temp.score));
            %             compression_scores(:,nn) = temp.score;
            %             continue;


            %             % RNN
            %             cfg=[];
            %             cfg.method='deep'; cfg.trn = result_tst;
            %             cfg.classifier_type = {'1d_cnn'}; %{'lda','domain_align',size(result_tst.trial)};
            %             cfg.perf_measure = 'acc';
            %             cfg.tst = result_trn; cfg.tst.trialinfo=result_trn.trialinfo(:,1); cfg.tst.trial = single(cfg.tst.trial);
            %             cfg.trn = result_tst; cfg.trn.trialinfo=result_tst.trialinfo(:,1);
            %             cfg.folds=nan; cfg.do_parallel=0;
            %             cfg.trn.trial = single(cfg.trn.trial);
            %             fidelity_pt(nn,:) = lv_classify(cfg);
            %             continue;

            % compressing temporal info. (mean, all pts.,fidelity/former)
            temp2 = repmat(result_tst.trialinfo(:,1),1,size(result_tst.trial,3)); result_tst.trialinfo=[];
            result_tst.trialinfo(:,1)=reshape(temp2,[],1);
            temp = permute(result_tst.trial, [2 1 3]);
            result_tst.trial = (reshape(temp, size(temp,1),[]))';
            id = find(isnan(sum(result_tst.trial,2)));  keep_id = find(~isnan(sum(result_tst.trial,2))); temp=nan(size(result_tst.trial,1),1);
            result_tst.trial(id,:,:)=[]; result_tst.trialinfo(id,:)=[];


            %             if size(result_trn.trial,2)==1, continue; end
            %             % lda beamformer
            %             cfg=[]; cfg.trn=result_tst; cfg.tst=result_trn;
            %             result = lv_mclass_beamformer(cfg);
            %             result_tst.trial = result.trn.trial; result_trn.trial = result.tst.trial;

            %             % choosing the best time pts based on the distance from zero
            %             cfg2=[]; cfg2.data(1,:,:) = abs(result_tst.trial);
            %             cfg2.method = 'fidelity'; %cfg2.percentile=1;
            %             [ good_trials ] = lv_reduce_trials(cfg2);

            % erps of the source cond_time_nn
            %             erp_trn(:,:,nn) = [squeeze(mean(result_tst.trial(result_tst.trialinfo(:,1)==1,:,:),1))';squeeze(mean(result_tst.trial(result_tst.trialinfo(:,1)==2,:,:),1))'];
            %             erp(:,:,nn) = [squeeze(mean(result_trn.trial(result_trn.trialinfo(:,1)==1,:,:),1))';squeeze(mean(result_trn.trial(result_trn.trialinfo(:,1)==2,:,:),1))'];
            %             continue;


            %             result_tst = lv_event_timelocked(result_tst, [0 1.5], good_trials); temp = result_tst.trialinfo;

            % %             % recurrence plot
            % %             cfg=[]; cfg.method = 'recurrence plot';
            % %             cfg.approach = 'correlation'; cfg.data = result_tst;
            % %             recurrence_mx(nn,:,:) = lv_post_classification(cfg);
            % %             continue;

            % % %             % getting fidelity from sleep sleep classification and using that wih wake
            % % %             cfg=[];
            % % %             cfg.method = 'axtime';
            % % %             cfg.classifier_type = {'lda','domain_align',size(result_tst.trial)}; cfg.perf_measure='dval'; % the size for domain align is overriden later in lv_classify
            % % %             cfg.trn = result_tst; cfg.trn.trialinfo=result_tst.trialinfo(:,1);
            % % %             cfg.folds=0; cfg.do_parallel=1;
            % % %             cfg.trn.trial = single(cfg.trn.trial);
            % % %             classification_dvals = lv_classify(cfg); %gives error be careful reduce data and don't use parallel ! and maybe use talldouble !
            % % %             cfg2=[]; cfg2.data(1,:,:) = classification_dvals;
            % % %             cfg2.method = 'fidelity'; cfg2.percentile=1;
            % % %             [ good_trials ] = lv_reduce_trials(cfg2);
            % % %
            % % %             % sleep wake classification with max fidelity point
            %             if sum(sum(squeeze(good_trials(1,:,:))))<50, [axtimeAcc(nn,:),axtimeAcc_notcentered(nn,:)]=deal(nan(1,length(result_trn.time)));  continue; end
            %             result_tst_temp = lv_event_timelocked(result_tst, [0 0], squeeze(good_trials(1,:,:)) ); temp = result_tst.trialinfo;
            %             result_tst_temp.trial = squeeze(result_tst_temp.trial);


            %             result_tst.trial = mean(squeeze(mean(result_tst.trial,2)),2);
            %             result_tst.trial = mean(result_tst.trial,3);


            % %             % sleep sleep classification
            % %             cfg=[];
            % %             cfg.method = 'axtime'; cfg.trn = result_tst; cfg.perf_measure = 'acc';
            % %             cfg.classifier_type = {'lda'}; %{'lda','domain_align',size(result_tst.trial)};
            % %             cfg.folds=5; cfg.do_parallel=1; cfg.trn.trialinfo=result_tst.trialinfo(:,1); cfg.trn.trial = single(cfg.trn.trial);
            % %             fidelity_pt2(nn,:) = lv_classify(cfg);
            % %             continue;

            %             % sleep sleep one pt
            %             result_tst.trial = zscore(result_tst.trial,[],1);
            %             Mdl = fitcdiscr(result_tst.trial, result_tst.trialinfo(:,1));
            %             [outclass,posterior] = predict(Mdl,result_tst.trial);
            %             fidelity_pt(nn,:) = sum(outclass==result_tst.trialinfo(:,1))./length(outclass);
            %             % getting the correct pts
            %             temp(keep_id) = outclass==result_tst.trialinfo(:,1);
            %             %temp(keep_id) = max(posterior,[],2) .* temp(keep_id); % putting the certainty here
            %             temp(temp==0)=nan;% temp(temp==0)=1;
            %             temp(isnan(temp))=0;
            %             res = reshape(temp, sz(1),sz(3)); clear temp; sleep_res = res;
            %             res = [nan(size(res,1),nearest(tst.data.time,0)-1) res nan(size(res,1),100)]; % to match the -4 4 length of trials
            %             temp(:,1,:) = res; temp=repmat(temp,1,size(tst.data.trial,2),1);


            %[corr_trl{nn,1},corr_time] = find(res==1);
            %             majority=[];
            %             for i=1:size(ss,1), majority(i,1) = sum(ss(i,:)==1)>sum(ss(i,:)==0); end % correct by majority of time
            %             corr_majority{nn,1} = find(majority==1);
            %             continue;

            % % %             % PAC of SO and spindles keep_id
            % % % %             tst.data.trial = tst.data.trial(keep_id,:,:);
            % % % %             tst.data.trialinfo = tst.data.trialinfo(keep_id,:);
            % % %             idx = cell2mat(arrayfun(@(x) (nearest(tst.data.time,x)), hold_trialinfo(:,6),'Un',0)); % variable length nans
            % % %
            % % %             keepidx = single(ones(size(tst.data.trial))); time0=nearest(tst.data.time,0);
            % % %             for i=1:size(tst.data.trial,1)
            % % %                 keepidx(i,:,idx(i)+1:end) = nan; keepidx(i,:,1:time0-1) = nan;
            % % %             end
            % % %             keepidx = temp .* keepidx;
            % % %             cfg=[]; cfg.data = tst.data;
            % % %             cfg.method = 'PAC'; cfg.keepidx = keepidx; % this is for the parts out of edges and nans of jittering
            % % %             [ result ] = lv_reduce_trials(cfg);
            % % % %             phase(nn,:) = result.phase; uniformity(nn,1) = result.uniformity;
            % % %             uniformity(nn,1) = result.coupling_strength;
            % % %             phase(nn,1) = result.phase_of_coupling; coupling_p(nn,1)=result.pval;

            %             % Temporal compression analysis .. sleep sleep and sleep wake correct pts.
            %             sz=size(result_trn.trial);
            %             result_trn.trial = zscore(result_trn.trial,[],1);
            %             result_trn.trial = permute(result_trn.trial, [2 1 3]);
            %             result_trn.trial = (reshape(result_trn.trial,size(result_trn.trial,1),[]))';
            %             lbls = repmat( result_trn.trialinfo(:,1),sz(3),1 );
            %
            %             [outclass,posterior] = predict(Mdl,result_trn.trial);
            %             temp=outclass==lbls;
            %             temp = reshape(temp,sz(1),sz(3));
            %             wake_res = temp;
            %             % sleep
            %             cfg=[]; cfg.data=sleep_res; cfg.method='compression_correct_pts';
            %             ratios = 10:701; % trusted reactivation lengths (50ms and longer)
            %             cfg.recurrence=1; cfg.ratios=ratios; % recurrence and the ratios that are considered as reactivation
            %             [temp,temp2] = lv_post_classification(cfg);
            %             lens_sleep = [lens_sleep temp];
            %             recurrence_distrib = [recurrence_distrib temp2];
            %
            %             % wake
            %             cfg=[]; cfg.data=wake_res;
            %             cfg.recurrence=0; cfg.method='compression_correct_pts';
            %             temp = lv_post_classification(cfg);
            %             lens = [lens temp];
            %
            %
            %             continue;



            % %             % RNN
            % %                         cfg=[];
            % %                         cfg.method='deep'; cfg.trn = result_tst;
            % %                         cfg.classifier_type = {'RNN'}; %{'lda','domain_align',size(result_tst.trial)};
            % %                         cfg.perf_measure = 'acc';
            % %                         cfg.folds=5; cfg.do_parallel=1; cfg.trn.trialinfo=result_tst.trialinfo(:,1); cfg.trn.trial = single(cfg.trn.trial);
            % %                         fidelity_pt(nn,:) = lv_classify(cfg);

            % % %             result_trn.trial = zscore(result_trn.trial,[],1); result_tst.trial = zscore(result_tst.trial,[],1); % zscoring
            % % % %             cfg=[]; cfg.method = 'timextime';
            % % % %             cfg.classifier_type = {'lda'}; %{'lda','domain_align',size(result_trn.trial)}; cfg.perf_measure='auc';
            % % % %             cfg.tst = result_trn; cfg.tst.trialinfo=result_trn.trialinfo(:,1);
            % % % %             cfg.trn = result_tst; cfg.trn.trialinfo=result_tst.trialinfo(:,1); cfg.folds=nan; cfg.do_parallel=1;
            % % % %             cfg.trn.trial = single(cfg.trn.trial); cfg.tst.trial = single(cfg.tst.trial);
            % % % %             classification_result{nn,cond} = lv_classify(cfg); % acc at every trn time
            % % %             clear temp;
            % % %


            % sleep wake classification
            result_trn.trial = zscore(result_trn.trial,[],1); result_tst.trial = zscore(result_tst.trial,[],1); % zscoring
            % % %             temp = mean(result_tst.trial,3); result_tst.trial=temp; clear temp;
            % %             %             temp = mean(mean(result_trn.trial,3),2); result_trn.trial=temp; clear temp;

            cfg=[]; cfg.method = 'timextime';
            cfg.classifier_type = {'lda'}; %{'lda','domain_align',size(result_trn.trial)};
            cfg.perf_measure = 'acc';
            cfg.tst = result_trn; cfg.tst.trialinfo=result_trn.trialinfo(:,1); cfg.tst.trial = single(cfg.tst.trial);
            cfg.trn = result_tst; cfg.trn.trialinfo=result_tst.trialinfo(:,1);
            cfg.folds=nan; cfg.do_parallel=1;
            cfg.trn.trial = single(cfg.trn.trial);
            fidelity_pt(nn,:) = lv_classify(cfg); % acc at every trn time

            continue;
            % %             % mvpa
            % %             rng(1); % to reproduce
            % %             cfg = [];
            % %             cfg.classifier      = 'lda';
            % %             cfg.metric          = 'auc'; % 'dval' distance from boundary 23taked
            % %             cfg.preprocess      =   {'undersample'};% de hat3ml zscoring le kol fold kaman 3shan mayb2ash fe leakage
            % %             cfg.hyperparameter           = [];
            % %             cfg.hyperparameter.lambda    = 'auto';
            % %             classification_result{nn,cond} = mv_classify_timextime(cfg,single(result_tst.trial) , result_tst.trialinfo(:,1), single(result_trn.trial) , result_trn.trialinfo(:,1)); %


            %             continue;

            % %                         % taking the peak trn time pt from wake
            % %                         cfg=[];
            % %                         cfg.method = 'axtime';
            % %                         cfg.classifier_type = {'lda'}; cfg.perf_measure='auc';
            % %                         cfg.trn = result_trn; cfg.trn.trialinfo=result_trn.trialinfo(:,1);
            % %                         cfg.folds=5; cfg.do_parallel=1; cfg.trn.trial = single(cfg.trn.trial);
            % %                         axTime_auc = lv_classify(cfg); [val,id] = max(axTime_auc);
            % %                         cfg=[]; cfg.latency=result_trn.time(id); result_trn = ft_selectdata(cfg,result_trn);
            % %                         % time locking to the most certain time pt of tst set trials
            % %                         % here we train with wake and test with sleep
            % %                         result_trn.trial = zscore(result_trn.trial,[],1); result_tst.trial = zscore(result_tst.trial,[],1); % zscoring
            % %                         cfg=[];
            % %                         cfg.method = 'timextime';
            % %                         cfg.classifier_type = {'lda'}; cfg.perf_measure='dval';
            % %                         cfg.trn = result_trn; cfg.trn.trialinfo=result_trn.trialinfo(:,1);
            % %                         cfg.tst = result_tst; cfg.tst.trialinfo=result_tst.trialinfo(:,1); cfg.folds=nan; cfg.do_parallel=1;
            % %                         cfg.trn.trial = single(cfg.trn.trial); cfg.tst.trial = single(cfg.tst.trial);
            % %                         classification_dvals = lv_classify(cfg); %gives error be careful reduce data and don't use parallel ! and maybe use talldouble !
            % %                         cfg2=[]; cfg2.data = (classification_dvals);
            % %                         cfg2.method = 'fidelity'; cfg2.percentile=1;
            % %                         [ good_trials ] = lv_reduce_trials(cfg2);
            % %
            % %                         for i=1:size(good_trials,1) % every trn pt will have a tst set so we cut every tst set
            % %                             if sum(sum(squeeze(good_trials(i,:,:))))<50, axtimeAcc(nn,i)=0.5; continue; end
            % %                             result_tst_temp = lv_event_timelocked(result_tst, [0 0], squeeze(good_trials(i,:,:)) ); temp = result_tst.trialinfo;
            % %                             result_tst_temp.trial = squeeze(result_tst_temp.trial);
            % %                             cfg=[]; cfg.method = 'timextime';
            % %                             cfg.classifier_type = {'lda'}; cfg.perf_measure='auc';
            % %                             cfg.trn = result_trn; cfg.trn.trialinfo=result_trn.trialinfo(:,1);
            % %                             cfg.trn.trial = squeeze(cfg.trn.trial(:,:,i));
            % %                             cfg.tst = result_tst_temp; cfg.tst.trialinfo=result_tst_temp.trialinfo(:,1); cfg.folds=nan; cfg.do_parallel=1;
            % %                             cfg.trn.trial = single(cfg.trn.trial); cfg.tst.trial = single(cfg.tst.trial);
            % %                             axtimeAcc(nn,i) = lv_classify(cfg); % acc at every trn time
            % %
            % %                         end
            % %                         continue;
            % %
            % %                         if size(result_tst.trial,1)<20, continue; end

            %             cfg=[];cfg.channel = lower({'C6','C4','C2','C1','C3','C5', 'CP5', 'CP3', 'CP1', 'CP2', 'CP4', 'CP6'}); % aoi
            %             %cfg.trials = find(sum(good_trials,2) > 0);
            %             result_tst=ft_selectdata(cfg,result_tst); result_tst.trialinfo = temp;





            %             % classification
            %             cfg=[];
            %             cfg.method = 'axtime';
            %             cfg.classifier_type = {'lda','domain_align',size(result_trn.trial)};
            %             cfg.perf_measure='auc';
            %             cfg.trn.trial = single(result_tst.trial); cfg.trn.trialinfo=result_tst.trialinfo(:,1);
            %             cfg.folds=5; cfg.do_parallel=1; %cfg.every_channel=1;
            %             sleep_classification_result(nn,:)  = (squeeze(lv_classify(cfg)))'; % ch_time inside sbj_cond

            % % %             cfg=[];
            % % %             cfg.method = 'axtime';
            % % %             cfg.classifier_type = {'lda'}; %,'domain_align',size(result_tst.trial)};%{'lda'};
            % % %             cfg.perf_measure='auc';
            % % %             cfg.trn.trial = single(result_trn.trial); cfg.trn.trialinfo=result_trn.trialinfo(:,1);
            % % %             cfg.folds=5; cfg.do_parallel=1; %cfg.every_channel=1;
            % % %             %                         classification_result{nn,cond} = lv_classify(cfg);
            % % %             wake_classification_result(nn,:)  = (squeeze(lv_classify(cfg)))'; % ch_time inside sbj_cond
            % % %
            % z_scoring
            %result_trn.trial = zscore(result_trn.trial,[],1); result_tst.trial = zscore(result_tst.trial,[],1); % zscoring

            %             cfg=[];
            %             cfg.method = 'timextime';
            %             cfg.perf_measure='auc';
            %             cfg.classifier_type = {'lda'};% {'lda','domain_align',size(result_trn.trial)};
            %             cfg.tst = result_trn; cfg.tst.trialinfo=result_trn.trialinfo(:,1);
            %             cfg.trn = result_tst; cfg.trn.trialinfo=result_tst.trialinfo(:,1); cfg.folds=nan; cfg.do_parallel=1;
            %             cfg.trn.trial = single(cfg.trn.trial); cfg.tst.trial = single(cfg.tst.trial);
            %             classification_result{nn,cond} = lv_classify(cfg); %gives error be careful reduce data and don't use parallel ! and maybe use talldouble !
            % %             % Rimannian
            % %             temp=[];
            % %             for i=1:size(result_trn.trial,1), temp(i,:) = reshape(cov( squeeze(result_trn.trial(i,:,:))' ), 1,[]); end
            % %             result_trn.trial = temp; temp=[];
            % %             for i=1:size(result_tst.trial,1), temp(i,:) = reshape(cov( squeeze(result_tst.trial(i,:,:))' ), 1,[]); end
            % %             result_tst.trial = temp; temp=[];
            % %
            % %             cfg=[];
            % %             cfg.method = 'timextime';
            % %             cfg.classifier_type = {'Riemannian','distance'}; % classifier, riemannian method
            % %             %cfg.tst = result_trn; cfg.tst.trialinfo=result_trn.trialinfo(:,1);
            % %             cfg.trn = result_tst; cfg.trn.trialinfo=result_tst.trialinfo(:,1); cfg.folds=5; cfg.do_parallel=0;
            % %             sleep_classification_result(nn,:) = lv_classify(cfg);
            % % %             % RSA
            % % %             cfg=[]; rng(1)
            % % %             result_trn_c1 = result_trn.trial(result_trn.trialinfo(:,1)==1,:,:);
            % % %             result_trn_c2 = result_trn.trial(result_trn.trialinfo(:,1)==2,:,:);
            % % %
            % % %             result_tst_c1 = result_tst.trial(result_tst.trialinfo(:,1)==1,:,:);
            % % %             result_tst_c2 = result_tst.trial(result_tst.trialinfo(:,1)==2,:,:);
            % % %
            % % %             cfg.data1.trial =  result_trn_c1; cfg.data2.trial =  result_tst_c1;
            % % %             cfg.dim_of_interest1 = 2; cfg.dim_of_interest2 = 2;
            % % %             cfg.handle = 'corr'; cfg.arguments = {'type' , 'pearson'};
            % % %             result_within_category1  = lv_trlxtrl_similarity(cfg);
            % % %             cfg.data1.trial =  result_trn_c2; cfg.data2.trial =  result_tst_c2;
            % % %             result_within_category2  = lv_trlxtrl_similarity(cfg);
            % % %             result_within_category(nn,:,:) = (result_within_category1+result_within_category2) ./2;
            % % %
            % % %             cfg.data1.trial =  result_trn_c1; cfg.data2.trial =  result_tst_c2;
            % % %             result_between_category1 = lv_trlxtrl_similarity(cfg);
            % % %             cfg.data1.trial =  result_trn_c2; cfg.data2.trial =  result_tst_c1;
            % % %             result_between_category2 = lv_trlxtrl_similarity(cfg);
            % % %             result_between_category(nn,:,:) = (result_between_category1+result_between_category2) ./2;
            % %             % 2d smoothing
            % %             K = ones(10); % that's 50ms smoothing
            % %             K = K.*(1/numel(K));
            % %             classification_result = cellfun(@(x) (conv2(x,K,'same')), classification_result,'Un',0);

            % %                         % mvpa
            % %                         rng(1); % to reproduce
            % %                         cfg = [];
            % %                         cfg.classifier      = 'lda';
            % %                         cfg.metric          = 'auc'; % 'dval' distance from boundary 23taked
            % %                         cfg.preprocess      =   {'undersample'};% de hat3ml zscoring le kol fold kaman 3shan mayb2ash fe leakage
            % %                         %cfg.k = 5;  cfg.method='timextime';
            % %                         cfg.hyperparameter           = [];
            % %                         cfg.hyperparameter.lambda    = 'auto';
            % %                         classification_result{nn,cond} = mv_classify_timextime(cfg,single(result_tst.trial) , result_tst.trialinfo(:,1), single(result_trn.trial) , result_trn.trialinfo(:,1)); %

            %  single(result_tst.trial) , result_tst.trialinfo(:,1));
            %             axtimeAcc(nn,:) = classification_result{nn,cond};
            % %             % axspace classifier
            % %             temp(:,1,:) = mean(result_tst.trial,3); result_tst.trial=temp; clear temp;
            % %             cfg=[];
            % %             cfg.method = 'axtime';
            % %             cfg.classifier_type = {'lda'};
            % %             cfg.do_parallel = 1;
            % %             %result_tst.trial = permute(result_tst.trial, [1 3 2]); % permute space and time
            % %
            % %             cfg.trn = result_tst; cfg.trn.trialinfo=result_tst.trialinfo(:,1); cfg.folds=5;
            % %             classification_axspace(:,nn) = lv_classify(cfg); % 1d curve in space




            %             % cfg = [];
            %             % cfg.layout = lay;   % this is the layout structure that you created with ft_prepare_layout
            %             % ft_layoutplot(cfg);
            %             ss = zeros(1,15); ss(4:6)=1;
            %             ft_plot_topo(anne_layout.pos(:,1), anne_layout.pos(:,2), ss, 'mask', anne_layout.mask,'outline' ,anne_layout.outline); % 2rsem 3leha el zvalues!
            %

        end
    end
    %     save classification_result classification_result % classification on every channel
    lv_pretty_errorbar(result_trn.time(nearest(result_trn.time,0):nearest(result_trn.time,1.1))...
        , axtimeAcc, (axtimeAcc*0)+0.5, 1);

    % topo visualisation
    [topo, ~]=ft_plot_topo(lv_layout.pos(:,1), lv_layout.pos(:,2), median(eigVectors,2),...
        'mask', lv_layout.mask,'outline' ,lv_layout.outline);


    result_tst.time = result_tst.time{1,1}; result_trn.time = result_trn.time{1,1};
    lv_pretty_errorbar(result_tst.time(nearest(result_tst.time,0):nearest(result_tst.time,2.5))...
        , squeeze(erp_trn(1,:,:))', squeeze(erp_trn(2,:,:))', 1);
    lv_pretty_errorbar(result_trn.time(nearest(result_trn.time,0):nearest(result_trn.time,1.1))...
        , squeeze(erp(1,:,:))', squeeze(erp(2,:,:))', 1);

    lv_pretty_errorbar(fidelity_pt(:,1),.25);
    res = lv_pretty_errorbar(result_trn.time, fidelity_pt, (fidelity_pt*0)+0.25, 1);
    res = lv_pretty_errorbar(result_trn.time, fidelity_pt(ids,:), (fidelity_pt(ids,:)*0)+0.25, 1);
    signi_pts = 34:139;
    c4 = fidelity_pt;
    hold on, lv_pretty_errorbar(result_trn.time, c4, (c4*0)+0.5, 99);

    lv_pretty_errorbar(1:58, spindle_likeli, (spindle_likeli*0)+0.5, 0);
    xticklabels(tst.data.label)
    [val,id]=max(mean(fidelity_pt, 1));
    lv_pretty_errorbar(fidelity_pt(:,id),0.25);

    for i=1:size(spindle_likeli,2), lv_vec2corr(fidelity_pt(:,id),spindle_likeli(:,i),'acc',tst.data.label(i)); end
    % for RSA
    classification_result=[];
    for i=1:size(result_within_category,1), classification_result{i,1}=squeeze(result_within_category(i,:,:)); end
    for i=1:size(result_between_category,1), classification_result{i,2}=squeeze(result_between_category(i,:,:)); end
    vschance=0;

    circ_plot(phase(uniformity<0.05),'pretty','bo',true,'linewidth',2,'color','b');  %phase plot of spindle max PAC
    [pval z] = circ_rtest(phase(uniformity<0.05)); hold on, title(['pval:' num2str(pval)]);
    [rho pval] = circ_corrcl(phase, max(fidelity_pt,[],2) ); % circular-linar correlation

    circ_plot(phase(coupling_p<0.05),'pretty','bo',true,'linewidth',2,'color','b');
    % phase analysis up vs down
    lv_pretty_errorbar(phase(:,1),.5);
    % canolty zscore of coupling strength to p-value
    p_one = normcdf(uniformity);  % because it looks at z<0 probabilities we take 1-pvalues to get the Ps for the positive strengths
    ppnt=find(uniformity>0&(1-p_one)<0.05);

    for i=1:size(wake_classification_result,1)
        figure, plot(result_trn.time, wake_classification_result(i,:)); title(['ppnt: ' num2str(sbj(i))]);
        if max(wake_classification_result(i,:))<0.85, disp(i);  end
    end
    classification_result([6 16 20 26 34 38],:)=[];

    lv_pretty_errorbar(result_trn.time, wake_classification_result, (wake_classification_result*0)+0.5, 1);
    lv_pretty_errorbar(result_tst.time, sleep_classification_result, (sleep_classification_result*0)+0.5, 1);

    axt_classification_result = cell2mat(cellfun(@(x) (mean(x,2)'),classification_result,'Un',0));% txt to axtime
    lv_pretty_errorbar(result_tst.time, axt_classification_result, (axt_classification_result*0)+0.5, 1);

    temp_classification_result = (cellfun(@(x) ((x(90:110,:))),classification_result,'Un',0)); %reducing time

    imagesc(result_tst.time,result_tst.time,squeeze(mean(recurrence_mx,1))); set(gca,'YDir','normal'); colorbar; caxis([-1 1])

    classification_result = classification_result(cell2mat(cellfun(@(x) (~isempty(x)),classification_result,'Un',0)));

    result_exp = classification_result(:,1);

    if vschance==1, result_chance = classification_result(:,1); else, result_chance = classification_result(:,2); end % chance or the control condition

    dat=[]; [cond1,cond2]=deal([]);
    for i=1:size(result_exp,1), cond1(i,:,:)=result_exp{i,1}; cond2(i,:,:)=result_chance{i,1}; end
    id = mod(1: size(cond1,1)+size(cond2,1), 2);
    if strcmp(cfg.method,'timextime')==1
        if vschance==1, dat.trial(id==0,1,:,:) = (cond2.*0)+0.25; else, dat.trial(id==0,1,:,:) = cond2; end % change to dat.trial(id==0,1,:,:) if many channels and txt
        dat.trial(id==1,1,:,:) = cond1;
        dat.time=result_trn.time; % dat.time=result_tst.time(nearest(result_tst.time,0):nearest(result_tst.time,1));
        dat.label={'z-stat'}; % freq should be trn
        dat.freq= result_tst.time( nearest(result_tst.time,0):nearest(result_tst.time,2.5));
        [ TF_struct ] = lv_tf(dat,1,1); %data, do_stats, do_plot
    else % axtime
        if vschance==1, dat.trial(id==0,:,:) = (cond2.*0)+0.50; else, dat.trial(id==0,:,:) = cond2; end
        dat.trial(id==1,:,:) = cond1;
        dat.time=result_trn.time; dat.label = result_trn.label;
        lv_erp(dat, 0, 1); %data, do_stats, do_plot
    end
    if ~isfield(cfg,'every_channel') && strcmp(cfg.method,'axtime')==1 % axtime stats when we use all channels and axtime
        cond1=squeeze(cond1);
        lv_pretty_errorbar(result_trn.time,cond1,(cond1.*0)+0.5, 2); end

    % the gcf has the result now save it for every run,, la de figure 100 w
    % kda el ids w b3d kda close el figures ba3d ma tesave ..
    if ~isempty(run_id), run_name = ['excel_' num2str(rund_id(run))]; else, run_name=num2str(run); end % if read from excel use the id to name the figures
    path = strcat(pathname, '\', ['run ' run_name]);
    fig_code = 100;
    figure(fig_code), saveas(gcf, strjoin([path ',Z' '.emf'])); saveas(gcf, strjoin([path ',Z' '.fig']));
    figure(fig_code+1), saveas(gcf, strjoin([path ',cond1' '.emf'])); saveas(gcf, strjoin([path ',cond1' '.fig']));
    figure(fig_code+2), saveas(gcf, strjoin([path ',cond2' '.emf'])); saveas(gcf, strjoin([path ',cond2' '.fig']));

    fprintf(['\n Finished run: ' run_name ' \n']);
end

% across space average plot as topo
[topo, ~]=ft_plot_topo(lv_layout.pos(:,1), lv_layout.pos(:,2), median(classification_axspace,2),...
    'mask', lv_layout.mask,'outline' ,lv_layout.outline);

figure; ft_plot_layout(lv_layout);

h = gca; h.XTickLabel = result_trn.label; xticks([1:length(result_trn.label)]); h.XTickLabelRotation = 90;
h.YTickLabel = result_trn.label; yticks([1:length(result_trn.label)]); h.YTickLabelRotation = 0;
h.TickLabelInterpreter = 'none';

% TF analysis
load tf_pw tf_pw
load TF_struct TF_struct
figure,
id1 = nearest(TF_struct.time,0):nearest(TF_struct.time,2.5); % not after 2.5 because of the jittering
TF_struct.time = TF_struct.time(id1);
id2 = nearest(TF_struct.freq,5):nearest(TF_struct.freq,30); 
TF_struct.freq = TF_struct.freq(id2);
temp = tf_pw(:,id2,id1);
b = imagesc(TF_struct.time, TF_struct.freq, squeeze(mean(temp,1)).*100 ); set(gca,'YDir','normal') 
xlabel('Time (sec.)', 'Interpreter','none');
ylabel('Frequency Hz', 'Interpreter','none');
h = colorbar; title(h,'Power');
caxis([0 25])

% ERP analysis
yyaxis right
id1 = nearest(tst.data.time,0):nearest(tst.data.time,2.5);
plot(tst.data.time(id1), mean(erp(:,id1),1), ...
    'black-')

% temporal compression
load lens lens
load lens_sleep lens_sleep
[vals_sleep,~] = histc(lens_sleep, 1:701); %linspace(min(lens),max(lens),nbin)
m = mode(lens_sleep); % make sure that the mode is the same as max
[val,id] = max(vals_sleep)
[val,id] = sort(vals_sleep,'descend');
bar(vals_sleep); % then inspect the histogram and see the values that are maximum of histo
ratios = [17 ]; ratios = [min(ratios):701];

[vals_wake,~] = histc(lens, 1:221); % wake
m = mode(lens);
[val,id] = max(vals_wake)
bar(vals_wake);
% 95 percentile
id95 = find(cumsum(val./sum(val))>0.95);
id95(1)

% wilson's compression
load compression_scores_smooth compression_scores_smooth
load compressions_smooth compressions_smooth
load compressionRatio compressionRatio
compression_scores = compression_scores_smooth;
avg_scores = mean(compression_scores,2);
avg_scores = avg_scores./sum(avg_scores);
sum(avg_scores'.*temp.compressionRatio) % weighted mean .. but affected by the value of the compressionRatio..
weighted_mean = mean(avg_scores.'*(temp.compressionRatio)',2); % this is from matlab answers
bar(mean(compression_scores,2))
% 95 percentile
id95 = find(cumsum(avg_scores,1)>0.95);
id95(1)
% compression vs dilation generically
lv_pretty_errorbar(mean(compression_scores(1:9,:),1)', mean(compression_scores(11:end,:),1)')

[~,id] = max(compression_scores,[],1);
sum( compressionRatio(id) < 1)
var(compressionRatio(id))

% significant scaling factors
load compression_scores compression_scores
load compressionRatio compressionRatio
lv_pretty_errorbar(compressionRatio, compression_scores', (compression_scores'.*0)+0.25,0);


% recurrence
dist_range = [min(recurrence_distrib) max(recurrence_distrib)];
[vals_recurr,~] = histc(recurrence_distrib, [dist_range(1):dist_range(2)]);
bar([dist_range(1):dist_range(2)],vals_recurr);
[val,id] = max(vals_recurr)
% ignoring 0bin and looking at the reactivation probabilty for each count
if dist_range(1)==0, vals_recurr=vals_recurr(2:end); dist_range=[min(recurrence_distrib)+1 max(recurrence_distrib)]; end
vals_recurr_norm = vals_recurr./sum(vals_recurr);
bar([dist_range(1):dist_range(2)],vals_recurr_norm.*100);
id=find(cumsum(sort(vals_recurr_norm,'descend'))>0.95) %95%
id(1)
%% PAC across time
id = pac_info(:,3)<0.05;
alpha=circ_mean(pac_info(id,2));
circ_plot(pac_info(id,2),'pretty','bo',true,'linewidth',2,'color','b')
circ_plot(pac_info(id,2),'hist',[],20,true,true,'linewidth',2,'color','b')
%% cross_freq coupling correlation with acc.
% pac_info(nn,ch,:) = [result.coupling_strength result.phase_of_coupling result.pval];
flags = squeeze(pac_info(:,:,3)); temp=flags;
flags(temp<0.05)=1; flags(temp>0.05)=0;
signi_ppnts = sum(flags,1);
acc = mean(fidelity_pt(:,34:139),2);

pac_signi = bsxfun(@times,pac_info,flags);
for i=1:size(flags,2)
    temp = flags(:,i);
    corr_stat = lv_vec2corr(pac_signi(temp~=0,i,1), acc(temp~=0),tst.data.label(i),'acc'); % correlation with coupling strength
    %     [corr_stat.R corr_stat.P]= circ_corrcc(pac_signi(temp~=0,i,2), acc(temp~=0)); % circular correlation with coupling phase
    if corr_stat.P>0.05, close(gcf); end, if corr_stat.P<0.05, disp(tst.data.label(i)); end
end
alpha=circ_mean( pac_signi(temp~=0,41,2) )
%% zvalue plot for more info about the difference

zval = zval_2d(cond1- ((cond1.*0)+0.5));


%% post-classification analyses

% correlation between classification and behaviour
% inside jittered_exp we have results.xlsx these before MRI and results_MRI_final.xlsx those are the MRI ppnts
% the labels are common and myDat are coming from the excel sheets and saved as mfiles
session = 3;
var = ['myDat_s' num2str(session)];
load (var);
load behav_lbls behav_lbls;
% temp = mean( fidelity_pt,1); [~,ii]=max(temp);
load fidelity_pt fidelity_pt;
pk = fidelity_pt(:,50); % 50 is the max point .. fidelity_pt(:,ii)
myDat = eval(var);
id = find(ismember(myDat(:,1), sbj(1:21))); % ids of sbj in excel
myDat = myDat(id,3:end); behav_lbls = behav_lbls(3:end);

% repeating for mri and then aggregate mydat
var = ['myDat_s' num2str(session) 'mri'];
load (var);
myDat_mri = eval(var);
id = find(ismember(myDat_mri(:,1), sbj(22:end))); % ids of sbj in excel
myDat_mri = myDat_mri(id,3:end);

myDat = [myDat ; myDat_mri];

for i=1:size(myDat,2)
    v = myDat(:,i); pk_temp=pk;
    if sum(isnan(v))>0 && sum(isnan(v))<length(v), pk_temp(isnan(v))=[]; v(isnan(v))=[]; end
    [Rt,Pt] = corr(pk_temp , v,'type','spearman'); % kendall sha3'al bardo
    if Pt<0.05
        stats = mo_vec2corr(pk_temp , v,'AUC',behav_lbls(i))
        title(['ppnts:' num2str(length(pk_temp))]);
    end
end
% for specific measures
measure = split('R_sss_late NR_sss_late');
myDat = myDat(:,ismember(behav_lbls,measure)); % then run the loop again
myDat(isnan(myDat(:,1)),:)=[];
lv_pretty_errorbar(myDat(:,1), myDat(:,2));

% MRI measures
load mri_measures mri_measures;
load mri_lbls mri_lbls;
id = ismember(mri_measures(:,1),sbj_MRI);
sbj_MRI_id = ismember(sbj_MRI,mri_measures(id,1)); sbj_MRI=sbj_MRI(sbj_MRI_id); % because not all ppnts got measures and not all of them in the excel

record=[]; pk = pk(22:end);
for i=1:length(sbj_MRI), id = ismember(mri_measures(:,1),sbj_MRI(i));
    if sum(id)>0, record=[record ; mri_measures(id,:) pk(i)]; end
end

myDat = record(:,2:end-1);
pk = record(:,end);
for i=1:size(myDat,2)
    v = myDat(:,i); pk_temp=pk;
    [Rt,Pt] = corr(pk_temp , v,'type','spearman');
    if Pt<0.05
        stats = mo_vec2corr(pk_temp , v,'AUC',mri_lbls(i))
    end
end

figure, imagesc(result.hfreq,result.lfreq,result.coupling_strength); set(gca,'YDir','normal'); colorbar; title('coupling strength (tort''s)')
figure, imagesc(result.hfreq,result.lfreq,result.pval<0.05); set(gca,'YDir','normal'); colorbar; title('significant pts')
%% behavioural improvement from TMR ?
% analysis of behvioural improvement across sessions caused by TMR
% preMRI
[re,nre,re_random,nre_random] = extract_blocks('myDat_s',sbj,2:4); % takes name of the dataset and sessions to analyse and returns the blocks aggregated
% MRI
[re2,nre2,re_random2,nre_random2] = extract_blocks('mri_myDat_s',sbj,2:4);
re=[re;re2]; nre=[nre;nre2]; re_random=[re_random;re_random2]; nre_random=[nre_random;nre_random2];

% for i=1:size(re,1), q1=prctile(re(i,:),5); re(i, re(i,:)>q1)=nan; end
% for i=1:size(nre,1), q1=prctile(nre(i,:),5); nre(i, nre(i,:)>q1)=nan; end
for i=1:size(re,1), q1=sort(re(i,:),'ascend'); re(i, ~ismember(re(i,:),q1(1:3)))=nan; end
for i=1:size(nre,1), q1=sort(nre(i,:),'ascend'); nre(i, ~ismember(nre(i,:),q1(1:3)))=nan; end

re = bsxfun(@minus,nanmedian(re_random,2),re); % vs. random as percentage change
nre = bsxfun(@minus,nanmedian(nre_random,2),nre);
lv_pretty_errorbar(nanmedian([re;re2],2),nanmedian([nre;nre2],2));

re1=re; nre1=nre; % post-sleep sesssions

% pre-sleep session1
[re,nre,re_random,nre_random] = extract_blocks('myDat_s',sbj,1);
[re2,nre2,re_random2,nre_random2] = extract_blocks('mri_myDat_s',sbj,1);
re=[re;re2]; nre=[nre;nre2]; re_random=[re_random;re_random2]; nre_random=[nre_random;nre_random2];
% for i=1:size(re,1), q1=prctile(re(i,:),5); re(i, re(i,:)>q1)=nan; end
% for i=1:size(nre,1), q1=prctile(nre(i,:),5); nre(i, nre(i,:)>q1)=nan; end
for i=1:size(re,1), q1=sort(re(i,:),'ascend'); re(i, ~ismember(re(i,:),q1(1:3)))=nan; end
for i=1:size(nre,1), q1=sort(nre(i,:),'ascend'); nre(i, ~ismember(nre(i,:),q1(1:3)))=nan; end

re = bsxfun(@minus,nanmedian(re_random,2),re); % vs. random as percentage change
nre = bsxfun(@minus,nanmedian(nre_random,2),nre);

% taking out session1 effect .. re1 should be sessions 2:4
ss = nanmedian([nanmedian(re,2)],2);
ids = find(ss<nanmedian(ss));

re1 = bsxfun(@minus,re1(ids,:),nanmedian(re(ids,:),2)); % gap post - gap pre
nre1 = bsxfun(@minus,nre1(ids,:),nanmedian(nre(ids,:),2));

lv_pretty_errorbar(nanmedian(re1,2),nanmedian(nre1,2));

% lv_pretty_errorbar(nanmedian(re(ids,:),2),nanmedian(nre(ids,:),2));

% lv_vec2corr(nanmedian(re1,2),acc,'R','acc')

% cfg=[]; cfg.method='pre-classification';
% cfg.data.trial=re'; temp = (lv_reduce_trials(cfg))'; temp(temp==0)=nan; re=re.*temp;
% cfg.data.trial=nre'; temp = (lv_reduce_trials(cfg))'; temp(temp==0)=nan; nre=nre.*temp;
% cfg.data.trial=re2'; temp = (lv_reduce_trials(cfg))'; temp(temp==0)=nan; re2=re2.*temp;
% cfg.data.trial=nre2'; temp = (lv_reduce_trials(cfg))'; temp(temp==0)=nan; nre2=nre2.*temp;
%% checking if presleep predicts postsleep performance
% regressing out presleep prediction of ovrenight improvement
re1=[re1;re21]; nre1=[nre1;nre21];
re=[re;re2]; nre=[nre;nre2];
mx = [nanmedian(re,2) nanmedian(re,2)-nanmedian(re1,2)];% ppntsx2
% mx=zscore(mx,[],1);
lv_vec2corr(mx(:,1),mx(:,2),'pre-sleep','change'); % will regress out the effect as in page 294 of Notes 9
% equation
x = mx(:,1);
y = (0.8071*x)-182.4;
y = mx(:,2)-y;
lv_vec2corr(x,y,'pre-sleep','Correctedchange');

R = y;
lv_pretty_errorbar(R,NR);

%% helping function
function data_parts = lv_check_parts(sbj ,type)
addDirs
if strcmp(type,'img')==1, data_parts=1; return; end
type = [mri_append 'part' num2str(sbj) '_' type '?'];
partIdx = strfind(type,'?');

for p=1:9 % 9 parts as maximum.. loop on this assumed parts until you know how many are there..
    type2 = type;
    type2(partIdx) = num2str(p);
    data_raw_file = (fullfile(rawdir, [type2 '.eeg']));
    if ~isfile(data_raw_file)
        data_parts = p-1;
        break;
    end

end

end

function data = lv_event_timelocked(data, duration, good_trials)
% to time lock the trials to an event with duration marking the length of
% the new trials and the starting point (because if the event happend just after TMR you maynot expect it to carry reactivation)
% after the beginning of the event .. the new trials and
% their new labels are returned in data
% good_trials trl_time and marks the event with ones
warning('assumed 200HZ sampling rate');
duration_samples = duration*200;

good_trials(:,1:duration_samples(1))=0;
good_trials(:,end-(duration_samples(2)-1):end)=0;

[events_row, events_col] = find(good_trials==1);

for i=1:length(events_row)
    new_event(i,:,:) = data.trial(events_row(i),:, events_col(i):events_col(i)+duration_samples(2) );
    new_label(i,:) = data.trialinfo(events_row(i),:);
end

data.trial = new_event; data.trialinfo = new_label; data.cfg=[];
class1_freq=sum(data.trialinfo(:,1)==1), class2_freq=sum(data.trialinfo(:,1)==2)
end

function result = lv_pac(cfg)

end

function zval = zval_2d(cond_difference)
% calculating the sample wise z-values using the difference between
% conditions given as: sbj_freq_time and returns the 2d z_vals
x = cond_difference;  zval=[]; % sbj_freq_time
for k=1:size(cond_difference,2) %freq
    for u=1:size(cond_difference,3) %time
        %         if sum(isnan(squeeze(x(:,k,u))))==0
        [~,~,stats] = signrank(squeeze(x(:,k,u)),0, 'method','approximate'); % you can change the stat function to something else here..(ex: ttest)
        zval(k,u) = stats.zval; % freq_time
        %         else zval(k,u)=nan; end
    end
end

end

function result = lv_mclass_beamformer(cfg)
% binary or multiclass beamforming
trn=cfg.trn;
tst=cfg.tst;

if size(trn.trial,2)==1, error('Data only got one feature, cannot do beamforming'); end

classes = unique(trn.trialinfo(:,1));
if numel(classes)==2, iterations=1; else iterations=numel(classes); end
for i=1:iterations
    % lda beamformer
    erp_pattern = mean( squeeze(mean(trn.trial(trn.trialinfo(:,1)==classes(i),:,:),1))...
        -squeeze(mean(trn.trial(trn.trialinfo(:,1)~=classes(i),:,:),1)), 2); % difference of classes erps is considered instead of cov_between .. and compressed in time
    epoched_data= permute(trn.trial,[2 3 1]); % 3D matrix (channels x time samples x epochs)
    [w1,temp(:,1,:),C1] = LDAbeamformer(erp_pattern,epoched_data);
    result.trn.trial(:,i,:) = temp; clear temp;

    % project the filter
    temp2 = permute(tst.trial,[2 3 1]);
    temp2 = w1' * reshape(temp2,size(temp2,1),[]);
    temp2 = reshape(temp2,size(tst.trial,3),size(tst.trial,1))';
    result.tst.trial(:,i,:) = temp2; clear temp2;
end

end


function [re,nre,re_random,nre_random] = extract_blocks(ds_name,sbj,sessions)
% takes name of the dataset and returns the blocks of all sessions
[re,nre,re_random,nre_random] = deal([]);
for i=1:24, R_pre{i,1} = ['R_pre_' num2str(i)]; NR_pre{i,1} = ['NR_pre_' num2str(i)]; end % blocks names
for i=1:2, R_random{i,1} = ['R_random_pre_' num2str(i)]; NR_random{i,1} = ['NR_random_pre_' num2str(i)]; end
if strcmp(ds_name,'myDat_s')==1, sbj=sbj(1:21); else, sbj=sbj(22:end); end
for session=sessions
    load behav_lbls behav_lbls;
    var = [ds_name num2str(session)];
    load (var); myDat = eval(var);
    id = find(ismember(myDat(:,1), sbj)); % ids of sbj in excel
    dat = myDat(id,3:end); behav_lbls = behav_lbls(3:end);

    re = [re dat(:, ismember(behav_lbls,R_pre))];% reactivated seq. blocks
    nre = [nre dat(:, ismember(behav_lbls,NR_pre))];
    re_random = [re_random dat(:, ismember(behav_lbls,R_random))];
    nre_random = [nre_random dat(:, ismember(behav_lbls,NR_random))];
end
end
%% this is a helping part that gets the previous/next trl's label and compares it to the current trl's label
% for nn=1:numel(sbj)
%     data_parts = lv_check_parts(sbj(nn),type);
%     cleaned_data = cell(1,data_parts);
%     for part=1:data_parts
%         cleaned_data{1,part} = lv_segment_filter_raw(sbj(nn),type, part, sleep_stage);
%     end
%     prev_lbl = cellfun(@(x) (x{1, 1}.prev), cleaned_data,'Un',0); prev_lbl = cell2mat(prev_lbl(:));
%     nxt_lbl = cellfun(@(x) (x{1, 1}.nxt), cleaned_data,'Un',0); nxt_lbl = cell2mat(nxt_lbl(:));
%     data = lv_load([cleaned_path num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)],'trial');
%     [~,~,rows] = intersect(data.trialinfo(:,[2 5]),prev_lbl(:,[2 5]), 'rows','stable'); prev_lbls = prev_lbl(rows,1);
%     [~,~,rows] = intersect(data.trialinfo(:,[2 5]),nxt_lbl(:,[2 5]), 'rows','stable'); nxt_lbls = nxt_lbl(rows,1);
%
%     [data_prev,data_nxt]=deal(data.trialinfo(:,1));
%     data_prev(prev_lbls==0,:)=[]; prev_lbls(prev_lbls==0)=[];
%     data_nxt(nxt_lbls==0,:)=[]; nxt_lbls(nxt_lbls==0)=[];
%
% %     data_prev(1)=[]; prev_lbls(1)=[];
% %     data_nxt()=[]; nxt_lbls()=[];
%
%     classes=unique(prev_lbls), if length(classes)~=4, error('lv: found classes are not 4 !'); end% the first two are aggregated together to be left hand and then the second two as right hand
%     prev_lbls( ismember(prev_lbls,classes(1:2)) ,1)  = 1; prev_lbls( ismember(prev_lbls,classes(3:4)) ,1)  = 2;
%     classes=unique(nxt_lbls), if length(classes)~=4, error('lv: found classes are not 4 !'); end% the first two are aggregated together to be left hand and then the second two as right hand
%     nxt_lbls( ismember(nxt_lbls,classes(1:2)) ,1)  = 1; nxt_lbls( ismember(nxt_lbls,classes(3:4)) ,1)  = 2;
%
%     classes=unique(data_prev), if length(classes)~=4, error('lv: found classes are not 4 !'); end% the first two are aggregated together to be left hand and then the second two as right hand
%     data_prev( ismember(data_prev,classes(1:2)) ,1)  = 1; data_prev( ismember(data_prev,classes(3:4)) ,1)  = 2;
%     classes=unique(data_nxt), if length(classes)~=4, error('lv: found classes are not 4 !'); end% the first two are aggregated together to be left hand and then the second two as right hand
%     data_nxt( ismember(data_nxt,classes(1:2)) ,1)  = 1; data_nxt( ismember(data_nxt,classes(3:4)) ,1)  = 2;
%
%     same_diff_prev(nn,:) = [mean(data_prev==prev_lbls) mean(data_prev~=prev_lbls)];
%     same_diff_nxt(nn,:) = [mean(data_nxt==nxt_lbls) mean(data_nxt~=nxt_lbls)];
%       'trial');
% end
%    lv_pretty_errorbar(same_diff_prev(:,1),same_diff_prev(:,2))
%    lv_pretty_errorbar(same_diff_nxt(:,1),same_diff_nxt(:,2))




%% TF analysis
function TFdat = do_tf(dat, baseline, frequencies) % takes 3d in .trial (trls_ch_time) and returns (ch_freq_time)

cfg              = [];
cfg.output       = 'pow';

cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:0.5:30; %linspace(frequencies(1),frequencies(end),2*(1+frequencies(end)-frequencies(1)));
cfg.t_ftimwin    = 5./cfg.foi;  % 5 cycles as a minimum to describe the frequency well
cfg.toi          = dat.time; % .time for max resolution .. to jump: window(1):0.1:window(2) this is just for visual smoothing
cfg.pad          ='nextpow2'; % rounds the maximum trial length up to the next power of 2
% cfg.keeptrials = 'yes';
TFdat = ft_freqanalysis(cfg, dat);


if ~isempty(baseline) && baseline(1)~=0
    cfg              = [];
    cfg.baseline     = [baseline(1) baseline(2)];
    cfg.baselinetype = 'relchange';
    [TFdat] = ft_freqbaseline(cfg, TFdat); % ch x freq x time
end

TFdat.trial = TFdat.powspctrm;

end




