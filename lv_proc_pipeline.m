%% main pipeline function that calls other functions
% add field trip path to load its default functions

% params
rawdir = 'D:\dataHub\lv_cleaning_erp_tf\all_data\Participant files\ppnt13';
addpath(['D:\sul''s', ' code\Toolboxes\fieldtrip-20190419'])
ft_defaults

type = 'sleep'; % 'sleep' 'img'
auto_cleaned_dir = [rawdir '\cleaned\'];

% naming convention: dataset_convention = 'part*_sleep?'; hypno_convention = 'psgHypno-part*_sleep?'; * for participant no. and ? for part no.
sleep_stage = 3;
stats_window = [-0.5 2.5]; % this is the time window for which we do the stats. for determining outliers and interpolation
fstruct = dir([rawdir '/part*.eeg']); h = struct2cell(fstruct);
sbj = unique(str2double( (cellfun(@(x) (strtok(x,['part_'])), (h(1,:))','Un', 0))' ));
pre_stim = 0.5;
post_stim = 2.5;
ref_ch = {''}; % {'TP9', 'TP10'};
required_sampling_rate = 200;


%% segmenting and cleaning
for nn=1:numel(sbj)
    data_parts = lv_check_parts(sbj(nn),type, rawdir);
    cleaned_data = cell(1,data_parts);
    for part=1:data_parts
        cleaned_data{1,part} = lv_segment_filter_raw(sbj(nn),type, part, sleep_stage, ref_ch,rawdir, pre_stim, post_stim, required_sampling_rate);
    end

    cleaned_data = cleaned_data(cell2mat(cellfun(@(x) (~isempty(x)),cleaned_data,'Un',0)));
    cleaned_data = ft_appenddata([],  cleaned_data{:}); % appending all parts together

    cleaned_data = lv_clean_segmented(cleaned_data, stats_window, sbj(nn));
    % convert to h5 for faster execution
    % lv_save([rawdir 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], cleaned_data, 'trial')
    cleaned_folder = fullfile(rawdir, 'cleaned');
    if ~exist(cleaned_folder, 'dir')
        mkdir(cleaned_folder);
    end
    new_path = fullfile(cleaned_folder, ['part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)]);
    lv_save(new_path, cleaned_data, 'trial');
end

%% manual artifact rejection for trials and channels .. giving trial numbers to be rejected then rejecting
% if we have different data and all is h5 we will need to use lv_save/load
% instead of the save/load

for nn=1:numel(sbj)
    fprintf(['\n Manual artifact rejection for trials and channels, subject: ' num2str(sbj(nn)) '\n']);
    data = lv_load([auto_cleaned_dir 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], 'trial');
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
    save([auto_cleaned_dir 'part' num2str(sbj(nn)) '_' type '_possible_bad_ch_N' num2str(sleep_stage)], 'possible_bad_ch', '-v7.3');
    save([auto_cleaned_dir 'part' num2str(sbj(nn)) '_' type '_manual_bad_trls_N' num2str(sleep_stage)], 'unique_id_bad_trls', '-v7.3');

    % rejecting trials
    mix_bad = unique_id_bad_trls;
    remaining_trls(nn,1)=size(data.trial,1)-size(mix_bad,1);
    [~,bad_trls] = intersect([data.trialinfo(:,2) data.trialinfo(:,5)] , mix_bad, 'rows'); % idx in trialinfo
    idx = 1:size(data.trial,1); idx(bad_trls)=[]; cfg = []; cfg.trials = idx;
    data = ft_selectdata(cfg, data);
    fprintf(['\n remaining_trls, subject: ' num2str(sbj(nn)) ', ' num2str(remaining_trls(nn,1)) '\n']);
    % saving
    cleaned_folder = fullfile(auto_cleaned_dir, 'final_cleaned_after_inspection');
    if ~exist(cleaned_folder, 'dir')
        mkdir(cleaned_folder);
    end
    lv_save([auto_cleaned_dir '\final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)], data, 'trial'); % cleaned data for further analyses
    fprintf(['\n Done. \n']);
end


%% rejecting arti. trials from data .. if more than one person was marking artifacts
% for nn=1:numel(sbj)
%     fprintf(['\n Manual artifact rejection MIXING, subject: ' num2str(sbj(nn)) '\n']);
%     %     load([preprocdir mri_append 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], 'cleaned_data'); % old save/load left commented because sleep data was saved this way
%     %     data = cleaned_data; clear cleaned_data; % data is a variable passed to the lv_manual_cleaning script not function to debug the code manually because it is manual inspection
%     data = lv_load([preprocdir mri_append 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], 'trial');
% 
%     % fixing IDs .. to make them unique because they repeat because of
%     % having different parts .. will encode it as real the trls and imaginary the part
%     if size(data.trialinfo,2)==4
%         part_shift = [0 find(diff(data.trialinfo(:,2))<0) length(data.trialinfo(:,2))]; no_prts=1:length(part_shift)-1;
%         if lv_check_parts(sbj(nn),type)~=length(no_prts), error('lv: parts mismatch between what is in data and what is on hard..'); end
%         for j=1:length(no_prts), temp_idx=part_shift(j)+1:part_shift(j+1); % to start from idx 1
%             data.trialinfo(temp_idx,5)= j; end % col. 5 for part no.
%     end
%     if sleep_stage==0 % wake
%         load([preprocdir 'manual_inspection_bad_trials\' mri_append 'part' num2str(sbj(nn)) '_' type '_manual_bad_trls_N' num2str(sleep_stage)], 'unique_id_bad_trls');
%         mix_bad = unique_id_bad_trls;
%     else
%         % loading martyna's bad trials and saved bad trials
%         %load([preprocdir 'manual_inspection_bad_trials\Martyna\data_cleaned\part' num2str(sbj(nn)) '_' type '_manual_bad_trls2_N' num2str(sleep_stage)], 'unique_id_bad_trls2');
%         load([preprocdir 'manual_inspection_bad_trials\' mri_append 'part' num2str(sbj(nn)) '_' type '_manual_bad_trls_N' num2str(sleep_stage)], 'unique_id_bad_trls');
%         mix_bad = unique_id_bad_trls; %union(unique_id_bad_trls2,unique_id_bad_trls,'rows'); % intersect/union applies unique automatically
%     end
%     remaining_trls(nn,1)=size(data.trial,1)-size(mix_bad,1);
%     [~,bad_trls] = intersect([data.trialinfo(:,2) data.trialinfo(:,5)] , mix_bad, 'rows'); % idx in trialinfo
% 
% 
%     idx = 1:size(data.trial,1); idx(bad_trls)=[]; cfg = []; cfg.trials = idx;
%     data = ft_selectdata(cfg, data);
% 
%     fprintf(['\n remaining_trls, subject: ' num2str(sbj(nn)) ', ' num2str(remaining_trls(nn,1)) '\n']);
% 
%     lv_save([preprocdir 'final_cleaned_after_inspection\' mri_append 'part' num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)], data, 'trial'); % cleaned data for further analyses
%     fprintf(['\n Done. \n']);
% end
% remaining_trls


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
conditions = 1; % if one condition then would repeat the data because there won't be stats.
erps_temp=[];
for nn=1:numel(sbj)
    fprintf(['\n ERP analysis, subject: ' num2str(sbj(nn)) '\n']);
    cleaned_data = lv_load([auto_cleaned_dir '\final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)],'trial');

    cfg=[]; cfg.latency=[-0.5 1.5]; cleaned_data=ft_selectdata(cfg,cleaned_data);% reducing to time of trial
    % putting trials of left hand together and for right hand as well
    if conditions == 2
        classes = unique(cleaned_data.trialinfo(:,1)); cleaned_data.trialinfo(ismember(cleaned_data.trialinfo(:,1),classes(1:2)),1 ) = 1; cleaned_data.trialinfo(ismember(cleaned_data.trialinfo(:,1),classes(3:4)),1 ) = 2;
    else
        cleaned_data.trial = [cleaned_data.trial ; cleaned_data.trial];
        cleaned_data.sampleinfo = [cleaned_data.sampleinfo ; cleaned_data.sampleinfo];
        cleaned_data.trialinfo = [cleaned_data.trialinfo ; cleaned_data.trialinfo];
        cleaned_data.trialinfo(1:size(cleaned_data.trialinfo,1)/2  ,1) = 1;
        cleaned_data.trialinfo((size(cleaned_data.trialinfo,1)/2)+1 : end) = 2; % ; (cleaned_data.trialinfo(:,1).*0)+2]
    end
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
    cleaned_data = lv_load([auto_cleaned_dir '\final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)],'trial');

    cfg=[]; cfg.latency=[-0.5 1.5]; cleaned_data=ft_selectdata(cfg,cleaned_data);% reducing to time of trial
    % putting trials of left hand together and for right hand as well
    if conditions == 2
        classes = unique(cleaned_data.trialinfo(:,1)); cleaned_data.trialinfo(ismember(cleaned_data.trialinfo(:,1),classes(1:2)),1 ) = 1; cleaned_data.trialinfo(ismember(cleaned_data.trialinfo(:,1),classes(3:4)),1 ) = 2;
    else
        cleaned_data.trial = [cleaned_data.trial ; cleaned_data.trial];
        cleaned_data.sampleinfo = [cleaned_data.sampleinfo ; cleaned_data.sampleinfo];
        cleaned_data.trialinfo = [cleaned_data.trialinfo ; cleaned_data.trialinfo];
        cleaned_data.trialinfo(1:size(cleaned_data.trialinfo,1)/2  ,1) = 1;
        cleaned_data.trialinfo((size(cleaned_data.trialinfo,1)/2)+1 : end) = 2; % ; (cleaned_data.trialinfo(:,1).*0)+2]
    end

    cleaned_data.baseline = [-0.5 0];
    [ TF_struct ] = lv_tf(cleaned_data, 0, 0); %data, do_stats, do_plot .. gets the TF in TF_struct.trial
    TF_temp = [TF_temp ; TF_struct.trial]; % erps_temp aggregates all the erps of different sbj
end

TF_temp_bl = TF_temp;  % with baseline
save TF_temp_bl TF_temp_bl;
save TF_temp TF_temp % 22sbj TF analyses

% group lvl TF
TF_struct.trial = TF_temp;
lv_tf(TF_struct, 0, 1); %data, do_stats, do_plot




%% helping function
function data_parts = lv_check_parts(sbj ,type, rawdir)
if strcmp(type,'img')==1, data_parts=1; return; end
type = ['part' num2str(sbj) '_' type '?'];
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
cfg.keeptrials = 'yes';
TFdat = ft_freqanalysis(cfg, dat);


if ~isempty(baseline) && baseline(1)~=0
    cfg              = [];
    cfg.baseline     = [baseline(1) baseline(2)];
    cfg.baselinetype = 'relchange';
    [TFdat] = ft_freqbaseline(cfg, TFdat); % ch x freq x time
end

TFdat.trial = TFdat.powspctrm;

end