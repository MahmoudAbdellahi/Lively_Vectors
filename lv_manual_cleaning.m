%% a script that takes data and performs different cleaning on it 
% it is a script to be able to change on the fly and run certain cleaning blocks

%% rejecting bad trials using stats and 1.5iqr and 25% channels to be sure of trials after interpolation 
cfg=[];
cfg.data.trial = mean(data.trial,3);
cfg.method = 'pre-classification';
[ good_trials1 ] = lv_reduce_trials(cfg);
cfg.data.trial = var(data.trial,0,3);
[ good_trials2 ] = lv_reduce_trials(cfg);
cfg.data.trial = max(data.trial,[],3);
[ good_trials3 ] = lv_reduce_trials(cfg);
cfg.data.trial = min(data.trial,[],3);
[ good_trials4 ] = lv_reduce_trials(cfg);
good_trials = good_trials1 & good_trials2 & good_trials3 & good_trials4;
possible_bad_trls = find( sum(~good_trials,2) > 0.25*size(data.trial,2)  ) 


%% rejecting bad trials without compressing time 
% in every channel the slice of trials_time then get the bad time points
% pts and then we sum them and see if their no. is higher than 1.5iqr+q3
% for >25% of channels we do not check low because now we sum the bad ones so the higher the worst and if no time pts are bad that is good!
cfg=[]; cfg.method = 'pre-classification'; bad_trls=[];
for i=1:size(data.trial,2)
    cfg.data.trial = squeeze(data.trial(:,i,:)); 
    [ good_trials_temp ] = lv_reduce_trials(cfg);
    trl_ch_bad = sum( ~logical(good_trials_temp) ,2);
    iqr_val = iqr(trl_ch_bad);  q3=prctile(trl_ch_bad,75);
    bad_trls(:,i)= trl_ch_bad > q3+(1.5*iqr_val);
end
possible_bad_trls_axtime = find(sum(bad_trls,2) > 0.25*size(data.trial,2))

%% identifying possible bad channel(s) with 1.5iqr (bad channel in majority of trials)
possible_bad_ch=[];
cfg = []; 
cfg.data.trial = mean(data.trial,3)'; % ch_trl, now for every trl get the channels representing outlires
cfg.method = 'pre-classification';
[ good_trials ] = lv_reduce_trials(cfg); ch  = (data.label( sum(good_trials,2) < round(size(data.trial,1)/2)  ));
fprintf(['\n Mean: possible bad channels are \n  ' cell2mat(ch(:)') '\n']);
possible_bad_ch{2,1} = ch;

cfg.data.trial = var(data.trial,0,3)'; 
[ good_trials ] = lv_reduce_trials(cfg); ch  = (data.label( sum(good_trials,2) < round(size(data.trial,1)/2)  ));
fprintf(['\n Variance: possible bad channels are \n  ' cell2mat(ch(:)') '\n']);
possible_bad_ch{2,2} = ch;

cfg.data.trial = max(data.trial,[],3)'; 
[ good_trials ] = lv_reduce_trials(cfg); ch  = (data.label( sum(good_trials,2) < round(size(data.trial,1)/2)  ));
fprintf(['\n Max: possible bad channels are \n  ' cell2mat(ch(:)') '\n']);
possible_bad_ch{2,3} = ch;

cfg.data.trial = min(data.trial,[],3)'; 
[ good_trials ] = lv_reduce_trials(cfg); ch  = (data.label( sum(good_trials,2) < round(size(data.trial,1)/2)  ));
fprintf(['\n Min: possible bad channels are \n  ' cell2mat(ch(:)') '\n']);
possible_bad_ch{2,4} = ch;

%% identifying possible bad channel with scattered channels (stats on time and median trials)
% this one compresses the channels so it is different from the one
% before because lv_reduce_trials will reduce the trials based on two
% thresholds but here we look at the general trend of abnormality based
% on avg. of trials .. so the one before based on no. trials but here based on strength of abnormality
% select scattered points then right click then create variable to
% store in variable: v1, v2, v3, v4
% difference between ft_ artifact rejection that here we can take the median
% mean
scat = median(mean(data.trial,3),1); % scattered channels
h = figure;
subplot(221), scatter(1:length(scat),scat); brush on, title('v1'), xlabel('channel'); ylabel('mean');
% var
scat = median(var(data.trial,0,3),1); % scattered channels
subplot(222), scatter(1:length(scat),scat); brush on, title('v2'), xlabel('channel'); ylabel('var');
% max
scat = median(max(data.trial,[],3),1); % scattered channels
subplot(223), scatter(1:length(scat),scat); brush on, title('v3'), xlabel('channel'); ylabel('max');
% min
scat = median(min(data.trial,[],3),1); % scattered channels
subplot(224), scatter(1:length(scat),scat); brush on, title('v4'), xlabel('channel'); ylabel('min');

% every measure got its own variable: v
% uiwait(h); % wait until figure is closed 
if exist('v1','var')==1, possible_bad_ch{2,5}= (data.label(v1(:,1))); clear v1; else possible_bad_ch{2,5}=[]; end
if exist('v2','var')==1, possible_bad_ch{2,6}= (data.label(v2(:,1))); clear v2; else possible_bad_ch{2,6}=[]; end
if exist('v3','var')==1, possible_bad_ch{2,7}= (data.label(v3(:,1))); clear v3; else possible_bad_ch{2,7}=[]; end
if exist('v4','var')==1, possible_bad_ch{2,8}= (data.label(v4(:,1))); clear v4; else possible_bad_ch{2,8}=[]; end

% storing in possible_bad_ch
possible_bad_ch(1,:) = {'mean','var','max','min','scattered_mean','scattered_var','scattered_max','scattered_min'};

%% suggesting bad trials using durations with amplitude threshold .. helps with inspection
% median and iqr based
temp=[];
for i=1:size(data.trial,1), temp{1,i} = squeeze(data.trial(i,:,:)); end % to 2d
dat = cell2mat(temp);

iquar = repmat(iqr(dat,2),1,size(dat,2));
q1=repmat(prctile(dat,25,2),1,size(dat,2)); q3=repmat(prctile(dat,75,2),1,size(dat,2));

bad_temp = (dat<q1-(2.*iquar) | dat>q3+(2.*iquar));
bad_segments = ( sum(bad_temp,1)>15 ); % 15 channels for more than 50ms .. change this with a value appropriate for your data
CC = bwconncomp(bad_segments); win_s = 100*0.05;
event = cell2mat(cellfun(@(x) (length(x)>win_s), CC.PixelIdxList, 'Un',0));
idx = CC.PixelIdxList(event)'; idx = cell2mat(idx);

trl_len = length(data.time);% samples to trials
iqr_bad_trls = unique(ceil(idx./trl_len))
    

possible_bad_trls_all = unique([possible_bad_trls ; possible_bad_trls_axtime ; iqr_bad_trls])
 
id_possible_bad_trls_all = [data.trialinfo(possible_bad_trls_all,2) data.trialinfo(possible_bad_trls_all,5)];
% trial browsing to see the rejected trials and then inspect without them.
temp_data = data; trl_len = length(temp_data.time);
begin = (1:trl_len:trl_len*size(temp_data.trial,1))'; temp_data.sampleinfo = [begin begin+trl_len-1];
cfg=[]; cfg.ylim=[-250 250]; % single trial based browsing
artf=ft_databrowser(cfg,temp_data);
bad_trl = ceil(artf.artfctdef.visual.artifact ./ length(temp_data.time)); % sample to trial
% 'bad_trl' will contain the indices of bad trials and you may want to
% aggregate them with 'possible_bad_trls_all' .. but here I chose to keep the
% trials from the analyses of 'possible_bad_trls_all' not the 'ft_databrowser'

idx = 1:size(data.trial,1); idx(possible_bad_trls_all)=[];
cfg = []; cfg.trials = idx;
iqr_data = ft_selectdata(cfg, data);

%% fieldtrip's manual visual inspection 
% ft_ data browser
% changing the original stored sampleinfo to be able to see the actual
% remaining trials not all the trials from the beginning of the pipeline..
temp_data = iqr_data;
trl_len = length(temp_data.time);
begin = (1:trl_len:trl_len*size(temp_data.trial,1))';
temp_data.sampleinfo = [begin begin+trl_len-1];

% cfg=[]; cfg.ylim=[-250 250]; % single trial based browsing
% artf=ft_databrowser(cfg,temp_data);
browser_cfg=[]; trls_to_show_together = 10; browser_cfg.ylim=[-250 250]; % continuous browsing
browser_cfg.blocksize = (abs(data.time(1))+abs(data.time(end))) * trls_to_show_together; % duration for cutting data
browser_cfg.continuous  = 'yes';
artf=ft_databrowser(browser_cfg,temp_data);

bad_trl = ceil(artf.artfctdef.visual.artifact ./ length(temp_data.time)); % ceil because we start from zero meaning that the first trial has samples<trial_len so the division will be zero


bad=[];
for i=1:size(bad_trl,1)
    bad = [bad bad_trl(i,1):bad_trl(i,2)];
end
bad = unique(bad);



unique_id_bad_trls = unique( [id_possible_bad_trls_all ; [temp_data.trialinfo(bad, 2) temp_data.trialinfo(bad, 5)] ] ,'rows'); % to save  





