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


%% parameters
rawdir = 'D:\dataHub\lv_cleaning_erp_tf\all_data\Participant files\ppnt13';
addpath('D:\codeHub\matlab\Toolboxes\fieldtrip-20190419') % fieldtrip-20190419 fieldtrip-20250107 fieldtrip-20250114
ft_defaults

type = 'sleep'; % 'sleep' 'img'
auto_cleaned_dir = [rawdir '\cleaned\'];

%% ERP and TF analyses
% the next two lines only for Marta's data because we don't have all the raw files here so
% would try to read the ppnts names from the cleaned
fstruct = dir([auto_cleaned_dir 'final_cleaned_after_inspection\part*.h5']); h = struct2cell(fstruct);
sbj = unique(str2double( (cellfun(@(x) (strtok(x,['part_'])), (h(1,:))','Un', 0))' ));

% would put that for one condition so won't be calling the function for two
% conditions
% parameters
channels = {'all'}; % all for average of all, or name of specific channel, or 
% if more than one channel it would be the average
tf_pw=[]; erp=[];
tf_baseline = [-0.4 0]; % -0.3 -0.1
tf_frequencies = [1 30];
erp_baseline = [-0.2 0];
temp = [];

for nn=1:numel(sbj)
    fprintf(['\n ERP analysis, subject: ' num2str(sbj(nn)) '\n']);
    cleaned_data = lv_load([auto_cleaned_dir '\final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_manual_cleaned_N' num2str(sleep_stage)],'trial');

    % time frequency analysis
    % choosing timepoint
    % cfg=[]; cfg.latency = [-0.4 2]; cleaned_data = ft_selectdata(cfg, cleaned_data); 
    % choosing channels
    cfg=[]; cfg.channel = channels; cleaned_data = ft_selectdata(cfg, cleaned_data);
    TF_struct= do_tf(cleaned_data, tf_baseline , tf_frequencies); % data, baseline, frequencies
    temp(1,:,:) = squeeze(mean(TF_struct.trial, 1)); % average of channels
    tf_pw(nn,:,:) = temp; % participant x frequency x time


    % ERP analysis
    baseline_id = nearest(cleaned_data.time,erp_baseline(1)):nearest(cleaned_data.time,erp_baseline(2)); % time
    baseline = (mean(cleaned_data.trial(:,:,baseline_id),3));
    baseline = repmat(baseline, 1,1,size(cleaned_data.trial,3));
    erp(nn,:) = mean(mean(cleaned_data.trial - baseline, 2),1);
end
save TF_struct TF_struct;
save tf_pw tf_pw;
save erp erp;
save cleaned_data cleaned_data;
%% plotting 
% tf analysis
load tf_pw tf_pw;
load TF_struct TF_struct;
load cleaned_data cleaned_data;

tf_frequencies = [5 30];
tf_duration = [-0.4 2];
erp_duration = [-0.4 2];

id1 = nearest(TF_struct.time,tf_duration(1)):nearest(TF_struct.time,tf_duration(2)); % not after 2.5 because of the jittering
TF_struct.time = TF_struct.time(id1);
id2 = nearest(TF_struct.freq, tf_frequencies(1)):nearest(TF_struct.freq, tf_frequencies(2)); 
TF_struct.freq = TF_struct.freq(id2);
temp = tf_pw(:,id2,id1);
b = imagesc(TF_struct.time, TF_struct.freq, squeeze(mean(temp,1)).*100 ); set(gca,'YDir','normal') 
xlabel('Time (sec.)', 'Interpreter','none');
ylabel('Frequency Hz', 'Interpreter','none');
h = colorbar; title(h,'Power');
caxis([-12 12])
% erp analysis
hold on,
load erp erp
yyaxis right
id1 = nearest(cleaned_data.time, erp_duration(1)):nearest(cleaned_data.time, erp_duration(2));
plot(cleaned_data.time(id1), mean(erp(:,id1),1), ...
    'black-')  
set(gca,'YColor','k'); % Change the right Axis's color to black 


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


% TF analysis
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
cfg.keeptrials   = 'no';
TFdat = ft_freqanalysis(cfg, dat);


if ~isempty(baseline) && baseline(1)~=0
    cfg              = [];
    cfg.baseline     = [baseline(1) baseline(2)];
    cfg.baselinetype = 'relchange';
    [TFdat] = ft_freqbaseline(cfg, TFdat); % ch x freq x time
end

% TFdat.powspctrm = mean(TFdat.powspctrm,1); % to get the average of trials and doing that here
% so that baseline is done on every trial

TFdat.trial = TFdat.powspctrm;

end