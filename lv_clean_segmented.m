function final_data_interp = lv_clean_segmented(segmented_data, stats_window, sbj)
% takes fieldtrip struct data and returns filtered fieldtrip struct 
% rejects outliers .. using min, max, var
% if one of the aoi (area of interest) is a bad channel
% then we reject that trial ... if >25% of channels are bad then reject the trial o.w. interpolate
% it using the neighboring channels.


if isfield(segmented_data,'lv_layout') % put the layout in segmented_data.lv_layout 
    lv_layout = segmented_data.lv_layout;
    aoi = lower(segmented_data.aoi);
else
    load lv_layout lv_layout; % defaults
    aoi = lower({'C6','C4','C2','C1','C3','C5', 'CP5', 'CP3', 'CP1', 'CP2', 'CP4', 'CP6'});
end
% area of interest (aoi)
bad_channels_threshold = 0.25;
aoi_idx = find(ismember( lower(segmented_data.label), aoi));

% to 3d if it isn't 3d .. 3d is the golden tensor
if length(size(segmented_data.trial)) < 3
    segmented_data.label = segmented_data.label;
    cfg             = [];
    cfg.keeptrials  = 'yes';
    data = ft_timelockanalysis(cfg, segmented_data);
else 
    data = segmented_data;
end

fprintf(['\n Rejecting outliers \n']);

% std_away = 2;
% Perc_threshold = 100*erf(std_away/sqrt(2)); % to keep the data around the mean approximately: 95% of the data with a room for error

cfg = []; 
cfg.latency = stats_window; % this is the time window for which we do the stats.
data_trialinfo = data.trialinfo;
[data] = ft_selectdata(cfg, data);
data.label = lower(data.label); segmented_data.label = lower(segmented_data.label);

%% getting the stats in time
all_var = squeeze(var(data.trial,0, 3)); % trl_ch
all_min = squeeze(min(data.trial,[], 3));
all_max = squeeze(max(data.trial,[], 3));
all_mean = squeeze(mean(data.trial, 3));

% calculating the thresholds and determining good trials
goodFlags_matrix = nan(size(all_var));
for i=1:size(all_var,2) 
    % interquartile range (IQR) based outliers rejection
    iqr_val = iqr(all_var(:,i)); q1=prctile(all_var(:,i),25); q3=prctile(all_var(:,i),75);
    var_threshold = [q1-(1.5*iqr_val)   q3+(1.5*iqr_val)]; % two sided threshold
    
    iqr_val = iqr(all_min(:,i)); q1=prctile(all_min(:,i),25); q3=prctile(all_min(:,i),75);
    min_threshold = [q1-(1.5*iqr_val)   q3+(1.5*iqr_val)];
    
    iqr_val = iqr(all_max(:,i)); q1=prctile(all_max(:,i),25); q3=prctile(all_max(:,i),75);
    max_threshold = [q1-(1.5*iqr_val)   q3+(1.5*iqr_val)];
    
    iqr_val = iqr(all_mean(:,i)); q1=prctile(all_mean(:,i),25); q3=prctile(all_mean(:,i),75);
    mean_threshold = [q1-(1.5*iqr_val)   q3+(1.5*iqr_val)];
    
    
    goodFlags_matrix(:,i) = all_var(:,i)>var_threshold(1) & all_var(:,i)<var_threshold(2) ...
        & all_min(:,i)>min_threshold(1) & all_min(:,i)<min_threshold(2) ...
        & all_max(:,i)>max_threshold(1) & all_max(:,i)<max_threshold(2) ...
        & all_mean(:,i)>mean_threshold(1) & all_mean(:,i)<mean_threshold(2);
end
% if we don't want to interpolate some channels that are very important to
% the analysis we reject the trial if bad on that important channel(s)
bad_aoi = sum(~goodFlags_matrix(:,aoi_idx),2) > length(aoi_idx)*(bad_channels_threshold);
goodFlags_matrix(bad_aoi , :) = 0;
fprintf(['\n ' num2str( sum(bad_aoi) ) ' trials rejected for being bad on aoi. \n']);

perfect_trls = find( sum(goodFlags_matrix,2)==size(goodFlags_matrix,2)); % all channels are good
all_good_trls = find( sum(goodFlags_matrix,2)./size(goodFlags_matrix,2) >= (1-bad_channels_threshold) );
bad_trls = 1:size(data.trial,1);
bad_trls(all_good_trls)=[];

fprintf(['\n ' num2str((length(all_good_trls)/size(all_var,1)) *100) '%% of trials accepted. \n']);

all_stats = [median(all_var,2) median(all_min,2) median(all_max,2)];


fprintf(['\n Interpolating ' num2str(length(all_good_trls)-length(perfect_trls)) ' trials. \n']);
% Find neighbours
cfg              = [];
cfg.method       = 'triangulation'; % the way it's going to connect channels indicating neighborhood
cfg.senstype     = 'EEG';
cfg.layout       = lv_layout;
cfg.feedback     = 'no';

neighbours	= ft_prepare_neighbours(cfg);

 
%% doing the interpolation by aggregating similar trials together
tic
data_interp = data;
[unique_trls,~,ic] = unique(goodFlags_matrix,'rows');
inter_co = 0; % interpolation count to make sure it was done correctly
for i=1:size(unique_trls,1)
    idx = find(ic == i);
    if any( ismember(union(perfect_trls,bad_trls), idx )  ), continue; end
    
    cfg                  = [];
    cfg.badchannel       = data.label( ~goodFlags_matrix(idx(1),:) ); 
    cfg.neighbours     	 = neighbours;
    cfg.layout           = lv_layout;
    cfg.method           = 'spline';
    cfg.feedback         = 'no';
    cfg.trackcallinfo	 = 'no';
    cfg.trials           = idx;
    
    interpAvg = ft_channelrepair(cfg ,data);
    
    % match the order of labels after the repair with the labels of data
    newIdx = lv_match_channels(data, interpAvg.label);
    
    if isfield(interpAvg,'avg'), data_interp.trial(idx,newIdx,:) = interpAvg.avg; inter_co=inter_co+length(idx); else
        data_interp.trial(idx,newIdx,:) = interpAvg.trial; inter_co=inter_co+length(idx); end
end
toc
cfg = [];
cfg.trials = all_good_trls;
[final_data_interp ] = ft_selectdata(cfg, data_interp);

final_data_interp.trialinfo = data_trialinfo(all_good_trls);

check_signals_plots(data_interp, perfect_trls, all_good_trls, all_stats, sbj)

if inter_co == length(all_good_trls)-length(perfect_trls) % to see if the interpolated trials are correctly considered
    fprintf(['\n Interpolation should be fine no. trials match interpolations. \n']);
else
    error('interpolated trials that should be interpoalted do not match what was actually interpolated!');
end
% % trial trial interpolation
% % data_interp = data;
% % sz = size(data.trial); % trls_ch_time
% % for i=1:sz(1)
% %     if ( sum(perfect_trls==i)+sum(bad_trls==i) )>0
% %         continue;
% %     end
% %
% %     % interpolate el data.trial
% %     cfg                  = [];
% %     cfg.badchannel       = data.label( ~goodFlags_matrix(i,:) );
% %     cfg.neighbourdist    = 4;
% %     cfg.neighbours     	 = neighbours;
% %     cfg.layout           = lv_layout;
% %     cfg.method           = 'spline';
% %     cfg.feedback         = 'no';
% %     cfg.trackcallinfo	 = 'no';
% %     cfg.trials           = i;
% %
% %     interpAvg = ft_channelrepair(cfg ,data);
% %
% %     % sanity check that the labels match
% %     lv_match_channels(data, interpAvg.label)
% %
% %     data_interp.trial(i,:,:) = interpAvg.avg;
% % end
% %
% %
% %
% % cfg = [];
% % cfg.trials = all_good_trls;
% % [final_data_interp] = ft_selectdata(cfg, data_interp);
% %
% % save final_data_interp final_data_interp

% checking the layout and plotting some data to see how it reflects on the
% topoplot, uses labels from the actual data and plots using the layout to make sure they match
% lv_match_channels(final_data_interp, lv_layout.label)
% datTemp = squeeze(final_data_interp.trial(1,:,1));
% datTemp(1) = 5000;
% figure,
% for i=1:length(final_data_interp.label)
% dat2 = circshift(datTemp, i-1);
% subplot(6,10,i),
% ft_plot_topo(lv_layout.pos(:,1), lv_layout.pos(:,2), dat2, 'mask', lv_layout.mask,'outline' ,lv_layout.outline);
% title(final_data_interp.label(i));
% end




end






function  check_signals_plots(data, perfect_trials, all_good_trls, all_stats, sbj)
% all channels some trials
figure('Name',['perfect trials, sbj:' num2str(sbj)]);
trls = perfect_trials( randi(length(perfect_trials),[1,10]) );
for i=1:length(trls)
    subplot(2,5,i);
    plot(data.time, squeeze(data.trial(trls(i),:,:)) );
    title({['var=' num2str(all_stats(trls(i),1))], ['min=' num2str(all_stats(trls(i),2))] ,...
        ['max=' num2str(all_stats(trls(i),3))]},'FontSize', 8);
end


figure('Name',['good trials, sbj:' num2str(sbj)]);
trls = all_good_trls( randi(length(all_good_trls),[1,10]) );
for i=1:length(trls)
    subplot(2,5,i);
    plot(data.time, squeeze(data.trial(trls(i),:,:)) );
    title({['var=' num2str(all_stats(trls(i),1))], ['min=' num2str(all_stats(trls(i),2))] ,...
        ['max=' num2str(all_stats(trls(i),3))]},'FontSize', 8);
end


all_bad_trls = 1:size(data.trial,1);
all_bad_trls(all_good_trls) = [];
figure('Name',['bad trials, sbj:' num2str(sbj)]);
trls = all_bad_trls( randi(length(all_bad_trls),[1,10]) );
for i=1:length(trls)
    subplot(2,5,i);
    plot(data.time, squeeze(data.trial(trls(i),:,:)) );
    title({['var=' num2str(all_stats(trls(i),1))], ['min=' num2str(all_stats(trls(i),2))] ,...
        ['max=' num2str(all_stats(trls(i),3))]},'FontSize', 8);
end

end



function newIdx = lv_match_channels(segmented_data, label)
% check that the labels are exactly the same

newIdx=[];
for i=1:length(label)
    newIdx = [newIdx ; find( ismember(segmented_data.label , label{i}) )];
end

idx = 1:length(label);
if any( idx' - newIdx  )
    warning('Channel mismatch');
end

end








