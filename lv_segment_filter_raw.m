function [ segmented_data ] = lv_segment_filter_raw(sbj,type, part, sleep_stage)
addDirs
fprintf(['\n Segmenting and filtering subject: ' num2str(sbj) '\n']);
typeEEG = type; 
if strcmp(type,'img')~=1 % imagery is in one part
    typeEEG = [mri_append 'part' num2str(sbj) '_' typeEEG '?'];
    typeEEG(strfind(typeEEG,'?')) = num2str(part);
else
    typeEEG = [mri_append 'part' num2str(sbj) '_' typeEEG]; rawdir = 'E:\MRI_jittered_exp\allIMG';% [datadir '/RawData/allImagery'];
end
data_raw    = (fullfile(rawdir, [typeEEG '.eeg']));
data_header = (fullfile(rawdir, [typeEEG '.vhdr']));

hdr = ft_read_header(data_raw);
nightinsamples = hdr.nSamples; % used in scoring
refIdx = find(ismember(hdr.label, [{'TP9'}, {'TP10'}])); % because ref. channels are always needed to get the value of any channel
refIdx = refIdx(:)';
if length(refIdx)~=2, error('lv: cannot find TP9 and TP10 !!'); end
b = waitbar(0, 'Segmenting and filtering data, please wait ...','WindowStyle','docked');
cfg             = [];
cfg.dataset     = data_raw;
cfg.channel = refIdx;
data_ref  = ft_preprocessing(cfg); % loading the reference once and append it to every channel later

skip_filtering = 0;
load_new_timelocking = 0;

for i=1: length(hdr.label)
    %% Filtering
    waitbar(i/length(hdr.label), b, ['Segmenting and filtering data, please wait ... channel: ' hdr.label(i)],'WindowStyle','docked');
    cfg             = [];
    cfg.dataset     = data_raw;
    cfg.channel = i;
    
    if any(refIdx==i)==1,  continue; end %ref
    
    tmp_data  = ft_preprocessing(cfg);
    tmp_data = ft_appenddata([], tmp_data,data_ref);
    
    if any(isnan(tmp_data.trial{1,1}(:))) % if nans in data then interpolate them
        waitfor(warndlg(['We have NANs in data in: ' num2str(sum(isnan(tmp_data.trial{1,1}(:)))) ' samples. press ok to interpolate with neighboring time points.']));
        tmp_data.trial{1,1} = fillmissing(tmp_data.trial{1,1} ,'linear',2);
    end
    
    cfg=[];
    cfg.reref       = 'yes';
    cfg.refchannel  = {'TP9', 'TP10'}; % the average of these two is used as the new reference, sometimes these become loose and that would be a problem
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [0.1 30]; % for new data sets consider FIR because it's better for offline analyses and not data driven as butter and IIR ! 
    % because IIR gets the filter kernel from data rather than sine so it can have the noise in data but that's ok
    cfg.bsfilter  = 'yes'; % band-stop method
    cfg.bsfreq    = [48 52];
    cfg.demean      = 'yes';
    
    if skip_filtering == 0
        tmp_data  = ft_preprocessing(cfg, tmp_data);
    end
    
    cfg2 = []; cfg2.channel = 1; % required channel is always the first one
    [data_proc{i}] = ft_selectdata(cfg2, tmp_data);
    clear tmp_data;
    
end
clear data_ref

data_proc = data_proc(~cellfun(@isempty, data_proc));
close(b);

for i=1:length(data_proc), data_proc{1,i}.label = lower(data_proc{1,i}.label); end % make the labels always lower case

%% Segmenting
if load_new_timelocking == 1 % if we use a different time locking and not TMR onset
    onsets = lv_load([num2str(sbj) '_time0']); % e.eg.[1000 4000 10200 25000]'; this mat should contain the vector of new onsets in samples
    [pre_stim,post_stim] = deal(4*data_proc{1,1}.fsample); % from seconds to samples
     
    trl=[];
    for i = 1:length(onsets)
        offset    = -pre_stim;  % number indicating the meaning of trlbegin
        trlbegin  = onsets(i) - pre_stim;
        trlend    = onsets(i) + post_stim;
        newtrl    = [trlbegin trlend offset];
        trl       = [trl; newtrl]; % store in the trl matrix
    end
    cfg     = [];
    cfg.trl = trl;
    data_def = ft_definetrial(cfg);
else
    [pre_stim,post_stim] = deal(4); % original
    cfg                     = [];
    cfg.dataset             = data_header;
    cfg.trialdef.eventtype  = 'Stimulus';
    cfg.trialdef.prestim    = pre_stim;
    cfg.trialdef.poststim   = post_stim;
    data_def = ft_definetrial(cfg);
end

for i=1:length(data_proc)
    segmented_data{i} = ft_redefinetrial(data_def, data_proc{i});
    data_proc(i)={1};
end
%% adding trial_lens field to get the original variable length of every trial.. important analysis after cleaning
% % % id0 = nearest(segmented_data{1, 1}.time{1, 1}, 0); % time 0 the actual onset
% % % for i=1:size(segmented_data{1, 1}.sampleinfo,1)
% % %     len = segmented_data{1, 1}.sampleinfo(i,1) : segmented_data{1, 1}.sampleinfo(i,2);
% % %     time0(i,1) = len(id0);
% % % end
% % % trial_lens = [diff(time0)./ segmented_data{1, 1}.fsample ; 3.5]; % the last 3.5 is for the last trial
% % % trial_lens = [segmented_data{1,1}.trialinfo (1:length(segmented_data{1,1}.trialinfo))' ...
% % %         segmented_data{1,1}.sampleinfo  repmat(part,length(segmented_data{1,1}.trialinfo),1) trial_lens];  
% % % % trial_lens(trial_lens(:,end)>3.5, end)=3.5; 
% % %  
% % % % exact segmenting by considering the jittering and not fixing trials' lengths ..
% % % % data_proc=hold_proc;
% % % cleaned_path = ['D:\sul''s code\Matt\sleep\erps\Organised\New exp\RawData\allSleep\data\cleaned\' ...
% % %                     'final_cleaned_after_inspection\MRI_part'];
% % % % cfg                     = [];
% % % % cfg.dataset             = data_header;
% % % % cfg.trialdef.eventtype  = 'Stimulus';
% % % % cfg.trialdef.prestim    = 0;
% % % % data_def = ft_definetrial(cfg);
% % % % data_def.trl(:,2)= [data_def.trl(2:end,1) ; data_def.trl(end,1)+(data_proc{1, 1}.fsample*3.5)]; % cutting trials based on jittering
% % % % 
% % % % pos_transition = find(data_def.trl(:,2)-data_def.trl(:,1)>data_proc{1, 1}.fsample*3.5); % long trials at the end of the sequence
% % % % data_def.trl(pos_transition,2) = data_def.trl(pos_transition,1)+data_proc{1, 1}.fsample*3.5;
% % % % 
% % % % bad_pos = find(data_def.trl(:,2)-data_def.trl(:,1)<data_proc{1, 1}.fsample*2.5); % very short trials that don't make sense
% % % % % data_def.trl(bad_pos,:) = [];
% % % % 
% % % % data_def.trl( data_def.trl(:,1)./data_proc{1,1}.fsample<4, :)=[]; 
% % % % if length(data_def.trl)~=length(segmented_data{1, 1}.trial)
% % % %     error('lv_ mismatch between segmented data triggers stored and jittered triggers');
% % % % end
% % % % trial_lens = (data_def.trl(:,2)-data_def.trl(:,1))./data_proc{1,1}.fsample;
% % % % trial_lens = [segmented_data{1,1}.trialinfo (1:length(segmented_data{1,1}.trialinfo))' ...
% % % %         segmented_data{1,1}.sampleinfo  repmat(part,length(segmented_data{1,1}.trialinfo),1) trial_lens2];  
% % % 
% % % load([cleaned_path num2str(sbj) '_' type '_manual_cleaned_N' num2str(sleep_stage)]);  
% % % 
% % % if ~isfield(other_data,'trial_lens') % then we entered here before so append 
% % %     other_data.trial_lens = trial_lens;  % then check the trials' lengths and remove the bad records
% % % else
% % %     other_data.trial_lens=[other_data.trial_lens;trial_lens];
% % % end
% % % 
% % % save([cleaned_path num2str(sbj) '_' type '_manual_cleaned_N' num2str(sleep_stage)], 'other_data', '-v7.3');  
% % % return;

% for i=1:length(data_proc)
%     segmented_data{i} = ft_redefinetrial(data_def, data_proc{i});
%     data_proc(i)={1};
% end


fprintf(['\n Finished segmenting and filtering subject: ' num2str(sbj) '\n']);

% unique trial identifier.. and the original sampleinfo, for SO and spindle marking later on
for i=1:length(segmented_data), segmented_data{1,i}.trialinfo = [segmented_data{1,i}.trialinfo (1:length(segmented_data{1,i}.trialinfo))' ...
        segmented_data{1,i}.sampleinfo  repmat(part,length(segmented_data{1,i}.trialinfo),1) ]; end  % last element is part no. to help uniquely identify the trial

% storing the original length of trials


% reverse looking at specific trials in here (original data) to see the label of their previous and next trials
% for i=1:length(segmented_data)
%     segmented_data{1,i}.prev = [0 ; segmented_data{1,i}.trialinfo(:,1)]; segmented_data{1,i}.prev(end,:)=[];
%     segmented_data{1,i}.prev = [segmented_data{1,i}.prev segmented_data{1,i}.trialinfo(:,2:end)];
%     segmented_data{1,i}.nxt = [segmented_data{1,i}.trialinfo(:,1) ; 0]; segmented_data{1,i}.nxt(1,:)=[]; % now the trials match look at the trials picked
%     segmented_data{1,i}.nxt = [segmented_data{1,i}.nxt segmented_data{1,i}.trialinfo(:,2:end)];
%     % and see their prev. and next
%
% end
% return

%% create layout
segmented_data = lv_create_layout(segmented_data, sbj);


%% sleep scoring
if sleep_stage~=0 % because 0 is imagery
    good_idx = lv_sleep_scoring_filtering(sbj, segmented_data, sleep_stage, type, part, nightinsamples);
else, good_idx = 1:length(segmented_data{1,1}.trial);
end
if isempty(good_idx), segmented_data=[]; return; end % should continue to the next part because this one is empty

cfg = []; cfg.trials = good_idx;
for i=1:length(segmented_data) % selecting trials in the correct sleep_stage and non-arousal
    segmented_data{1,i} = ft_selectdata(cfg,segmented_data{1,i});
end
% return;
% % % counting stims.
% % fprintf('----------------------------------------------------------------');
% % fprintf(['\n ppnt: ' num2str(sbj) ',part ' num2str(part)])
% % sz = size(segmented_data{1, 1}.trialinfo  ,1);
% % if segmented_data{1, 1}.trialinfo(2,1)-segmented_data{1, 1}.trialinfo(1,1)>10 || segmented_data{1, 1}.trialinfo(2,1)-segmented_data{1, 1}.trialinfo(1,1)==0, segmented_data=sz-1, else, segmented_data=sz, end
% % return
%% downsampling the data
cfgr = [];
fs = 200;
% will unify the sampling rate by specifying the time points that we want
% instead of the sampling rate: cfgr.resamplefs = 200; because when we have
% different sampling rates before doing this even we downsample it will give different no. samples
x = (-pre_stim:1/fs:post_stim); % even better to be like this because you will guarantee that the time isn't skewed it was (-4:3.995) but now it is symmetric
cfgr.time = repmat(mat2cell(x,1,length(x)), 1,length(segmented_data{1,1}.trial));

for i=1:length(segmented_data)
    segmented_data{1,i} = ft_resampledata(cfgr, segmented_data{1,i});
end
fprintf(['\n The sampling rate is now set to: ' num2str(segmented_data{1,1}.fsample) ' samples/sec. \n']);

%% appending different channels in one wholesome dataset
segmented_data = ft_appenddata([], segmented_data{:});
fprintf(['\n Data from different channels were appended successfully. \n']);

end









%% helping functions

function segmented_data = lv_create_layout(segmented_data, sbj)
% matches the labels between the actual data and the layout .. removes ref
% gnd emg eog. and keeps only those that exist in the default layout
%%

fprintf(['\n Creating layout. \n']);
labels =  ( cellfun(@(x) (cell2mat(x.label)),segmented_data,'Un',0) )';

cfg = [];
cfg.layout = 'easycapM1.mat';
lay = ft_prepare_layout(cfg);

lay.label = lower(lay.label); labels = lower(labels);  % because the labels were sometimes capital and sometimes small

select = ismember(lay.label,  labels);

select(ismember( lay.label ,lower({'TP9','TP10','CPz'}))) = 0; % removing the new references

lay.pos = lay.pos(select,:); lay.height = lay.height(select); lay.label = lay.label(select); lay.width = lay.width(select);


cfg = [];
cfg.layout = lay;   % this is the layout structure that you created with ft_prepare_layout
ft_layoutplot(cfg);
set(gcf,'Name',['layout sbj:' num2str(sbj)],'NumberTitle','off');

% if the current sbj has the same layout as the prev. one then it is safe to proceed
load lv_layout lv_layout;
for i=1:length(lay.label)
    if find(ismember(lv_layout.label,lay.label(i)))~=i, error('lv: layout of sbj does not match with the previous sbjs!'); end
end


lv_layout = lay;
save lv_layout lv_layout;


keep = ismember(labels,lv_layout.label);
segmented_data = segmented_data(1,keep);

segmented_data = lv_match_channels(segmented_data, lv_layout.label); % change channels order to match the layout

% checking the layout and plotting some data to see how it reflects on the topoplot
% load lv_layout lv_layout
%
% dat = [0 ones(1,length(lv_layout.label)-1)];
% figure,
% for i=1:length(lv_layout.label)
% dat2 = circshift(dat, i-1);
% subplot(6,10,i),
% ft_plot_topo(lv_layout.pos(:,1), lv_layout.pos(:,2), dat2, 'mask', lv_layout.mask,'outline' ,lv_layout.outline);
% title(lv_layout.label(i));
% end

fprintf(['\n Finished.. removed eog, emg, ref, and gnd to match with the default layout. \n']);
end


%%
function segmented_data2 = lv_match_channels(segmented_data, label)
% changes the order of the data inside the struct segmented_data
% such that it matches the channel order in label

labels =  ( cellfun(@(x) (cell2mat(x.label)),segmented_data,'Un',0) )';

newIdx=[];
for i=1:length(label), newIdx = [newIdx ; find( ismember( labels  , label{i} ) )]; end

segmented_data2 = segmented_data(1,newIdx);

end

%%
function good_idx = lv_sleep_scoring_filtering(sbj, segmented_data, stage, type, part, nightinsamples)
addDirs 
type = ['psgHypno-' mri_append 'part' num2str(sbj) '_' type '?'];
type(strfind(type,'?')) = num2str(part);

rawhypno = load(fullfile(rawdir, [type '.mat']));
if size(rawhypno.dat,1)>1, error('lv: scoring should be one vector !'); end
rawscore = rawhypno.dat;
for i=1:length(rawhypno.arousals), if isempty(rawhypno.arousals{i,1}), rawhypno.arousals(i,:)=[]; end, end
rawarousal = (cell2mat(rawhypno.arousals));
dat = segmented_data{1,1};

full_scored_samples = 10.*ones(1,  nightinsamples);

% rawscore 30 sec. scores
scored_samples=[]; arousal_samples=[]; good_idx = []; Epochtime = 30; arousal_idx=[];
for i=1:length(rawscore), scored_samples = [scored_samples repmat(rawscore(1,i),1, Epochtime*dat.fsample)]; end % the last unscored samples are not considered
full_scored_samples(1:length(scored_samples)) = scored_samples;
for i=1:size(rawarousal,1), arousal_samples = [arousal_samples rawarousal(i,1)*dat.fsample:rawarousal(i,2)*dat.fsample]; end

nonstage_samples = find(full_scored_samples  ~= stage);
remove_samples = unique(union(arousal_samples , nonstage_samples));

% just stage checker for the no. of trials at each stage to know how many times they woke up etc.
% adding arousals to the scoring
arousal_samples(round(arousal_samples)==0)=1; % because sometimes the first one is time0 so that's the first sample..
full_scored_samples(round(arousal_samples))=99;
stages = [0 1 2 3 5 10 99];
stage_trl_count = zeros(1,length(stages));


wake_samples = find(full_scored_samples==0); wake_counter=0; % just for wake
% dat.sampleinfo(end,:)=[];
agg_trls=[]; % hypno of trls
for i=1:size(dat.sampleinfo,1)
    if isempty( intersect(dat.sampleinfo(i,1):dat.sampleinfo(i,2) , remove_samples) )
        good_idx = [good_idx i];
    end
    
    if ~isempty( intersect(dat.sampleinfo(i,1):dat.sampleinfo(i,2) , arousal_samples) ) % just for the checking plot of arousals
        arousal_idx = [arousal_idx i];
        %continue; % because if the trl is an arousal we won't care what the 30sec. apoch says as it will be arousal trial
    end
    
    % just stage checker for the no. of trials at each stage to know how many times they woke up etc.
    temp = full_scored_samples(dat.sampleinfo(i,1):dat.sampleinfo(i,2));
    agg_trls = [agg_trls temp];
    stage_majority = mode(temp); % majority falling in which stage?
    stage_trl_count(stages==stage_majority) = ...
        stage_trl_count(stages==stage_majority)+1;
    
    if length(intersect(dat.sampleinfo(i,1):dat.sampleinfo(i,2),wake_samples)) > ...
            length(dat.sampleinfo(i,1):dat.sampleinfo(i,2))/2
        wake_counter = wake_counter+1;
    end
    %     for k=1:length(stage)
    %         if length( intersect(dat.sampleinfo(i,1):dat.sampleinfo(i,2) , stage{1,k}) )==length(dat.sampleinfo(i,1):dat.sampleinfo(i,2))
    %             stage_trl_idx{1,k} = [stage_trl_idx{1,k} i]; % purely in stage
    %             break;
    %         elseif ...
    %             ~isempty( intersect(dat.sampleinfo(i,1):dat.sampleinfo(i,2) , stage{1,k}) )&&(k==1 || k==6)
    %             bad_transition_trl_idx = bad_transition_trl_idx+1;% bad trial has movement-10/wake-0/arousal
    %             break;
    %         elseif length( intersect(dat.sampleinfo(i,1):dat.sampleinfo(i,2) , stage{1,k}) ) >  ...
    %             length(dat.sampleinfo(i,1):dat.sampleinfo(i,2))/2
    %             stage_trl_idx{1,k} = [stage_trl_idx{1,k} i]; % transition with non-arousal
    %             break;
    %         end
    %     end
    
end
plot(agg_trls)
lv_save(['D:\sul''s code\Matt\sleep\erps\Organised\New exp\stages\' mri_append 'sbj' num2str(sbj) '_part' num2str(part) ], stage_trl_count);

fprintf(['\n Finished scoring subject: ' num2str(sbj) ', no. trials: ' num2str(length(good_idx)) '\n']);

if isempty(good_idx), return; end % no trials in the needed sleep stage

% checking some good trials vs some bad trials
% some good trials
% % timeax = dat.time{1,1};
% % randtrls = randi(length(dat.trial),[1 10]);
% % figure,
% % for i=1:length(randtrls)
% %     subplot(2,5,i),
% %     plot(timeax, dat.trial{1,randtrls(i)} );
% %     title(['non-arousal trial:' num2str(randtrls(i)) ]);
% % end
% % % some arousals
% % randtrls = arousal_idx( randperm(length(arousal_idx)) );
% % if length(arousal_idx)>10, randtrls = randtrls(1:10); end
% % figure,
% % for i=1:length(randtrls)
% %     subplot(2,5,i),
% %     plot(timeax, dat.trial{1,randtrls(i)} );
% %     title(['arousal trial:' num2str(randtrls(i)) ]);
% % end



end



