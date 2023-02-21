function [ segmented_data ] = lv_segment_filter_raw_lateral(sbj,type, part, sleep_stage)
addDirs
fprintf(['\n Segmenting and filtering subject: ' num2str(sbj) '\n']);
typeEEG = type;
if strcmp(type,'img')~=1 % imagery is in one part
    typeEEG = ['part' num2str(sbj) '_' typeEEG '?'];
    typeEEG(strfind(typeEEG,'?')) = num2str(part);
else
    typeEEG = ['part' num2str(sbj) '_' typeEEG]; rawdir = [datadir '/RawData/allImagery'];
end
data_raw    = (fullfile(rawdir, [typeEEG '.eeg']));
data_header = (fullfile(rawdir, [typeEEG '.vhdr']));

hdr = ft_read_header(data_raw);
nightinsamples = hdr.nSamples; % used in scoring
% refIdx = find(ismember(hdr.label, [{'TP9'}, {'TP10'}])); % because ref. channels are always needed to get the value of any channel
% refIdx = refIdx(:)';
% if length(refIdx)~=2, error('lv: cannot find TP9 and TP10 !!'); end
b = waitbar(0, 'Segmenting and filtering data, please wait ...','WindowStyle','docked');
% cfg             = [];
% cfg.dataset     = data_raw;
% cfg.channel = refIdx;
% data_ref  = ft_preprocessing(cfg); % loading the reference once and append it to every channel later
 
for i=1:length(hdr.label)
    %% Filtering
    waitbar(i/length(hdr.label), b, ['Segmenting and filtering data, please wait ... channel: ' hdr.label(i)],'WindowStyle','docked');
    cfg             = [];
    cfg.dataset     = data_raw;
    cfg.channel = i; 
    
%     if any(refIdx==i)==1,  continue; end %ref
    
    tmp_data  = ft_preprocessing(cfg); 
%     tmp_data = ft_appenddata([], tmp_data,data_ref);
    
    if any(isnan(tmp_data.trial{1,1}(:))) % if nans in data then interpolate them
        waitfor(warndlg(['We have NANs in data in: ' num2str(sum(isnan(tmp_data.trial{1,1}(:)))) ' samples. press ok to interpolate with neighboring time points.'])); 
        tmp_data.trial{1,1} = fillmissing(tmp_data.trial{1,1} ,'linear',2); 
    end 

    cfg=[];
%     cfg.reref       = 'yes';
%     cfg.refchannel  = {'TP9', 'TP10'}; % the average of these two is used as the new reference, sometimes these become loose and that would be a problem
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [0.1 30]; % for new data sets consider FIR because it's better for offline analyses and not data driven as butter and IIR !
    cfg.bsfilter  = 'yes'; % band-stop method
    cfg.bsfreq    = [48 52];
    cfg.demean      = 'yes';
    
    tmp_data  = ft_preprocessing(cfg, tmp_data);
     
    cfg2 = []; cfg2.channel = 1; % required channel is always the first one
    [data_proc{i}] = ft_selectdata(cfg2, tmp_data);
    clear tmp_data;
    
end
clear data_ref

data_proc = data_proc(~cellfun(@isempty, data_proc));
close(b);
 
for i=1:length(data_proc), data_proc{1,i}.label = lower(data_proc{1,i}.label); end % make the labels always lower case

%% Segmenting
[pre_stim,post_stim] = deal(4);

cfg                     = [];
cfg.dataset             = data_header;
cfg.trialdef.eventtype  = 'Stimulus';
cfg.trialdef.prestim    = pre_stim;
cfg.trialdef.poststim   = post_stim;
data_def = ft_definetrial(cfg);

for i=1:length(data_proc)
    segmented_data{i} = ft_redefinetrial(data_def, data_proc{i});  
    data_proc(i)={1};
end


fprintf(['\n Finished segmenting and filtering subject: ' num2str(sbj) '\n']);

% unique trial identifier.. and the original sampleinfo, for SO and spindle marking later on
for i=1:length(segmented_data), segmented_data{1,i}.trialinfo = [segmented_data{1,i}.trialinfo (1:length(segmented_data{1,i}.trialinfo))' ...
        segmented_data{1,i}.sampleinfo  repmat(part,length(segmented_data{1,i}.trialinfo),1) ]; end  % last element is part no. to help uniquely identify the trial

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



% erp on raw test
%% timelockanalysis puts trials into a more convenient 3-d matrix
cfg             = [];
cfg.keeptrials  = 'yes';
im = ft_timelockanalysis(cfg, segmented_data);

clabel_hand_train = im.trialinfo(:,1);
clabel_hand_train( ismember(clabel_hand_train ,[1 2]) ) = 1;
clabel_hand_train( ismember(clabel_hand_train,[3 4]) ) = 2;

% cfg=[];
% cfg.reref       = 'yes';
% cfg.refchannel  = {'cz'};  
% im  = ft_preprocessing(cfg, im);

cfg=[];
% cfg.channel = {'c5','cp3'};
cfg.channel = {'c5','c6'};
im=ft_selectdata(cfg,im);

% erps
im.trial = im.trial(:,1,:) - im.trial(:,2,:);
cfg=[]; cfg.latency=[-1.1 2.2]; im = ft_selectdata(cfg,im);

erp1 = squeeze(mean(im.trial(clabel_hand_train==1,:,:),2));
erp2 = squeeze(mean(im.trial(clabel_hand_train==2,:,:),2)); 
errorbar(im.time, mean(erp1,1) , std(erp1,[],1)./sqrt(size(erp1,1))); hold on,
errorbar(im.time, mean(erp2,1) , std(erp2,[],1)./sqrt(size(erp2,1)));
 


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

type = ['psgHypno-part' num2str(sbj) '_' type '?'];
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

for i=1:size(dat.sampleinfo,1)
    if isempty( intersect(dat.sampleinfo(i,1):dat.sampleinfo(i,2) , remove_samples) )
        good_idx = [good_idx i];
    end
    
    if ~isempty( intersect(dat.sampleinfo(i,1):dat.sampleinfo(i,2) , arousal_samples) ) % just for the checking plot of arousals
        arousal_idx = [arousal_idx i];
    end
end


fprintf(['\n Finished scoring subject: ' num2str(sbj) ', no. trials: ' num2str(length(good_idx)) '\n']);

if isempty(good_idx), return; end % no trials in the neeeded sleep stage

% checking some good trials vs some bad trials
% some good trials
timeax = dat.time{1,1};
randtrls = randi(length(dat.trial),[1 10]);
figure,
for i=1:length(randtrls)
    subplot(2,5,i),
    plot(timeax, dat.trial{1,randtrls(i)} );
    title(['non-arousal trial:' num2str(randtrls(i)) ]);
end
% some arousals
randtrls = arousal_idx( randperm(length(arousal_idx)) );
if length(arousal_idx)>10, randtrls = randtrls(1:10); end
figure,
for i=1:length(randtrls)
    subplot(2,5,i),
    plot(timeax, dat.trial{1,randtrls(i)} );
    title(['arousal trial:' num2str(randtrls(i)) ]);
end



end



