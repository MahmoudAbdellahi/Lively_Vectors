function lv_plot_stats(data, conds, do_stats, pos )
% explores the dimensions of the data .. stats and plotting ..

% and does the stats if needed and plots
% the stats for lvl1(difference) or 2(corrected) and uses topographgy of data.lv_layout  
% if it is logical to plot on topo ...

% WARNING: when correcting we use the whole window as the clustering window so cut the data
% beforehand. in time, freq


elle m5alene msh 3ayz 2a3adel hwa en el curves w el 2d mo5talfeen ana msh damn awe awe eno nafs el concept ya3ne ene lw 3amalt el freq=1 kda b2a 1d
bs mashy 7ases eno kda bnsba kbera lakn still bardo lw kan kda msh hayb2a el donia 7lwa awe ene 25aleha 7aga wa7da le sabab en hayb2a fe ta3del fe el plotting
la2n el curves btb2a errors lel 2 curves lakn el figures 3ade kan kol figure lwa7daha ,,,
%% Curves
if length(size(data.trial))==3 %rpt_ch_time
    if isfield(data,'trialinfo')
        % sbj lvl
        for i=1:size(data.trial,2)
            axes('Position',[pos(i,1) pos(i,2)  0.05 0.05]) % awel 2 el x w el y aken el figure nafsaha mn 0 l7ad 1 w da el makan elle feh el figure
            % sbj lvl no correction:
            mo_pretty_errorbar(data.time, squeeze(data.trial(data.trialinfo==conds(1),i,:)),squeeze(data.trial(data.trialinfo==conds(2),i,:)), 99);
            title(data.label(i));
        end
    else
        % group lvl
        cluster_in_space = 0;
        parametric = 0; % change for different stat.
        cond1_idx = mod(1:size(data.trial,1),2);
        if cluster_in_space==0
            for i=1:size(data.trial,2)
                axes('Position',[pos(i,1) pos(i,2)  0.05 0.05]) % awel 2 el x w el y aken el figure nafsaha mn 0 l7ad 1 w da el makan elle feh el figure
                if do_stats==0
                    mo_pretty_errorbar(data.time,squeeze(data.trial(logical(cond1_idx) ,i,:)),squeeze(data.trial(logical(~cond1_idx) ,i,:)), 99);
                else
                    disp( string(data.label(i)) );
                    if parametric==1
                        mo_pretty_errorbar(data.time,squeeze(data.trial(logical(cond1_idx) ,i,:)),squeeze(data.trial(logical(~cond1_idx) ,i,:)), 1);
                    else
                        mo_pretty_errorbar(data.time,squeeze(data.trial(logical(cond1_idx) ,i,:)),squeeze(data.trial(logical(~cond1_idx) ,i,:)), 2);
                    end
                end
                title(data.label(i));
            end
            
        else
            % I am leaving correct_in_space however it will give an error
            % because it is not done in fieldtrip for ERPs I only find it
            % for TF analysis, TODO: correct in space by compressing time topoplot bs ya3ne: ft_topoplotER
            mo_pretty_errorbar(data.time,squeeze(data.trial(logical(cond1_idx) ,:,:)),squeeze(data.trial(logical(~cond1_idx) ,:,:)), 1, cluster_in_space);
        end
    end
    if size(data,2)>1 % many channels
        % interactive plotting, send the data in the figure handle and plot the
        % avg. of chosen channels
        load lv_layout lv_layout;
        figure
        cfg = [];
        cfg.layout  = lv_layout;
        lay = ft_prepare_layout(cfg);
        ft_plot_layout(lay)
        % add the required guidata
        info       = guidata(gcf);
        info.x     = lv_layout.pos(:,1);
        info.y     = lv_layout.pos(:,2);
        info.label = lv_layout.label;
        info.data = data;
        info.data.pos = pos;
        info.do_stats = do_stats;
        if ~isempty(parametric), info.parametric = parametric; end
        
        guidata(gcf, info)
        %   GUIDATA(H, DATA) stores the specified data in the figure's
        %   application data.  H is a handle that identifies the figure
        
        set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', @compress_channels, 'event', 'WindowButtonDownFcn'})
        set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', @compress_channels, 'event', 'WindowButtonDownFcn'})
        set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', @compress_channels, 'event', 'WindowButtonDownFcn'})
    end
end

%% TF like data
if length(size(data.trial)) > 3 %takes 'rpt_chan_freq_time'
    %takes 'rpt_chan_freq_time' in tfr.powspctrm
    data = rmfield(data ,'trial');
    
    if size(data.powspctrm,1)>2 % group lvl
        % doing the stats
        parametric = 0; % change for different stat.
        if parametric==1, warning('lv: performing parametric t-stat on sample lvl as stats');
            perform_TF_correction(data, pos);
        else, warning('lv: performing non-parametric wilcoxon on sample lvl as stats');
            perform_TF_correction_nonparametric(data, pos);
        end
    else
        % perform difference, no correction we just have one sbj
        perform_difference(data, data.lv_layout.pos);
    end
    
    if size(data,2)>1 % many channels
        % making it interactive because the averaging wasn't handled in fieldtrip in the same way as averaging and then correcting again which we want
        % interactive plotting, send the data in the figure handle and plot the avg. of chosen channels
        %                 layout_name = lv_tune_params('type the layout name to be loaded:', 'lv_layout');
        %                 tmp = load(layout_name);
        %                 lv_layout = eval(['tmp.' num2str(layout_name)]);
        lv_layout = data.lv_layout;
        
        figure
        cfg = [];
        cfg.layout  = lv_layout;
        lay = ft_prepare_layout(cfg);
        ft_plot_layout(lay)
        % add the required guidata
        info       = guidata(gcf);
        info.x     = lv_layout.pos(:,1);
        info.y     = lv_layout.pos(:,2);
        info.label = lv_layout.label;
        info.data = data;
        info.do_stats = do_stats;
        info.parametric = parametric;
        
        guidata(gcf, info)
        %   GUIDATA(H, DATA) stores the specified data in the figure's
        %   application data.  H is a handle that identifies the figure
        
        set(gcf, 'WindowButtonDownFcn',   {@ft_select_channel, 'multiple', true, 'callback', @compress_channels, 'event', 'WindowButtonDownFcn'})
        set(gcf, 'WindowButtonUpFcn',     {@ft_select_channel, 'multiple', true, 'callback', @compress_channels, 'event', 'WindowButtonDownFcn'})
        set(gcf, 'WindowButtonMotionFcn', {@ft_select_channel, 'multiple', true, 'callback', @compress_channels, 'event', 'WindowButtonDownFcn'})
    end
    
    
end

end




function   compress_channels(channels) % @ft_select_channel .. of fieldtrip will pass the channels to this handle here and the data is stored in the figure
% handle so that we can use it here.

info        = guidata(gcf); %   DATA = GUIDATA(H) returns previously stored data, or an empty

selected_channels = ismember(info.data.label,channels);

%% Curves
if isfield(info.data,'trial') %rpt_time
    compressed_temp = squeeze(mean(info.data.trial(:,selected_channels,:), 2)); % % trl_time
    
    data = info.data;
    data.trial = compressed_temp;
    figure,
    if isfield(data,'trialinfo')
        % sbj lvl
        conds = unique(data.trialinfo); if length(conds)>2, error('data has more than two conditions!!'); end
        % shouldn't make stats for sbj lvl
        mo_pretty_errorbar(data.time, squeeze(data.trial(data.trialinfo==conds(1),:)),squeeze(data.trial(data.trialinfo==conds(2),:)), 99);
        
        str = strjoin(string(data.label(selected_channels)'),', ');
        if sum(selected_channels)<length(data.label), title(strjoin(['mean(' str ')'])); else
            title('Grand average of all channels'); end
    else
        % group lvl
        cond1_idx = mod(1:size(data.trial,1),2);
        if info.do_stats==0
            mo_pretty_errorbar(data.time,squeeze(data.trial(logical(cond1_idx) ,:)),squeeze(data.trial(logical(~cond1_idx) ,:)), 99);
        else
            if info.parametric==1
                mo_pretty_errorbar(data.time,squeeze(data.trial(logical(cond1_idx) ,:)),squeeze(data.trial(logical(~cond1_idx) ,:)), 1);
            else
                mo_pretty_errorbar(data.time,squeeze(data.trial(logical(cond1_idx) ,:)),squeeze(data.trial(logical(~cond1_idx) ,:)), 2);
            end
        end
        str = strjoin(string(data.label(selected_channels)'),', ');
        if sum(selected_channels)<length(data.label), title(strjoin(['mean(' str ')'])); else
            title('Grand average of all channels'); end
    end
end



%% TF like data
if isfield(info.data,'powspctrm')
    data = info.data;
    data.powspctrm = mean(data.powspctrm(:,selected_channels,:,:),2);
    
    if length(data.label)==sum(selected_channels), data.label = {'Grand average of all channels'}; else
        str = strjoin(string(data.label(selected_channels)'),', ');
        data.label = { char(strjoin(['mean(' str ')'])) }; end
    if size(data.powspctrm,1)>2 % group lvl
        if info.parametric==1
            perform_TF_correction(data, []);
        else
            perform_TF_correction_nonparametric(data, []);
        end
    else
        % sbj lvl just plot the difference
        perform_difference(data, []);
    end
end



end



function  perform_TF_correction(tfr, pos)
% performs the cluster permutation test and uses plot2d for plotting in pos. of channels
cond1S = tfr; cond2S = tfr;

cond1_idx = mod(1:size(tfr.powspctrm,1),2);
cond1S.powspctrm = tfr.powspctrm(logical(cond1_idx), :,:,:); %'subj_chan_freq_time'
cond2S.powspctrm = tfr.powspctrm(~logical(cond1_idx), :,:,:);
cfg             = [];
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.clusteralpha = 0.05;        % alpha for thresholding the t statistic

cfg.latency     =  [tfr.time(1) tfr.time(end)];% expected time of effect assumed the whole time window
cfg.frequency   =  [tfr.freq(1) tfr.freq(end)];
% Setting up the cluster test
cfg.correctm            = 'cluster';
cfg.clusterstatistic    = 'maxsum';
cfg.minnbchan           = 0;
cfg.alpha               = 0.025; % since two-tailed
cfg.tail                = 0;
cfg.clustertail         = 0;
cfg.numrandomization    = 1000;

cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';

cluster_in_space = 1;
if cluster_in_space==0
    cfg.neighbours       = [];
else
    % if you want to cluster in space then look at the figure with channels don't choose some channels
    load lv_layout lv_layout; cfg_temp = [];
    cfg_temp.method   = 'triangulation'; cfg_temp.senstype = 'EEG'; cfg_temp.layout = lv_layout;
    cfg.neighbours = ft_prepare_neighbours(cfg_temp);
end
% Set up design matrix
% Row #1: subject nr (repeat for both conditions)
% Row #2: condition 1, then condition 2
sbj = 1:size(cond1S.powspctrm,1);
design = zeros(2,2*numel(sbj));
design(1,:)= [1:numel(sbj), 1:numel(sbj)];
design(2,1:numel(sbj))        = 1;
design(2,numel(sbj)+1:2*numel(sbj)) = 2;

cfg.design      = design;
cfg.uvar        = 1;  % unit variable in 1st row
cfg.ivar        = 2;  % independent variable in 2nd row (cond)


cond1S.cumtapcnt=[];
cond1S.grad=[];
cond1S.cfg=[];
cond1S.dimord = 'subj_chan_freq_time';

cond2S.cumtapcnt=[];
cond2S.grad=[];
cond2S.cfg=[];
cond2S.dimord = 'subj_chan_freq_time';

if cluster_in_space==0
    labels = cond1S.label;
    hold_cond1S_powspctrm = cond1S.powspctrm;
    hold_cond2S_powspctrm = cond2S.powspctrm;
    temp = cond1S; temp.powspctrm = []; % to hold the masked results for plotting
    for i=1:size(hold_cond1S_powspctrm,2) % on every channel
        cond1S.label = labels(i); cond2S.label = labels(i);
        cond1S.powspctrm = hold_cond1S_powspctrm(:,i,:,:);
        cond2S.powspctrm = hold_cond2S_powspctrm(:,i,:,:);
        
        [statsStruct] = ft_freqstatistics(cfg, cond1S, cond2S);
        
        % plot masked 2d on positions with condition1 because it is the one
        % of interest HENA!!!!!!!!!!!!!!!
        statsStruct.mask = double(statsStruct.mask);
        temp.powspctrm(i,:,:) = squeeze(statsStruct.stat); % to make temp ch_freq_time
        mask = squeeze(statsStruct.mask);
    end
    plot2d(temp ,'t-stat, corrected with cluser permutation',pos, mask )
else
    
    % WARNING: trace with actual data !
    cfg_main = cfg;
    [statsStruct] = ft_freqstatistics(cfg, cond1S, cond2S);
    % plotting significant clusters' tstat
    cfg = [];
    cfg.layout = lv_layout;
    cfg.parameter = 'powspctrm';
    cfg.maskparameter = 'mask';
    cfg.maskstyle = 'outline'; cfg.showlabels = 'yes';
    
    % getting the difference
    tfr_difference = tfr;
    tfr_difference.powspctrm = (cond1S.powspctrm-cond2S.powspctrm) ./ (cond1S.powspctrm+cond2S.powspctrm);
    tfr_difference.powspctrm = squeeze(mean(tfr_difference.powspctrm,1));% since many sbj we will get the avg. difference
    tfr_difference.mask = statsStruct.mask;
    
    ft_multiplotTFR(cfg, tfr_difference);
    c = colorbar('location', 'southoutside');
    c.Label.String = 'Power ratio (cond1 - cond2)/(cond1 + cond2)';
    
    
    figure, % topo of the effect, avg over time and freq. so you should limit the data to the intersting part before running it
    %     ft_topoplotER(cfg, tfr_difference);
    cfg_main.avgovertime = 'yes'; % if you avg like that then if the TF itself has nans you won't get significant clusters in space because you're averaging nans
    cfg_main.avgoverfreq = 'yes';
    [stats_topo] = ft_freqstatistics(cfg_main, cond1S, cond2S);
    figure,
    cfg = [];
    cfg.alpha  = 0.05;
    cfg.parameter = 'stat';
    cfg.layout    = lv_layout ;
    ft_clusterplot(cfg, stats_topo);
    title('compressed over all time and all freq. with alpha=0.05');
    %     ft_movieplotTFR(cfg, stats_topo);
end
% % check stats of single channel to see that actually when channels are put
% % together they have impact on each other !!!
% chech_ch = find(ismember(tfr.label,'Cz'));
% cond1S.powspctrm = cond1S.powspctrm(:,chech_ch,:,:);
% cond2S.powspctrm = cond2S.powspctrm(:,chech_ch,:,:);
% cond1S.label = {'mu'}; cond2S.label = {'mu'};
% [statsSingleChannel] = ft_freqstatistics(cfg, cond1S, cond2S );
% temp = statsStruct.mask(chech_ch,:,:);
% any(temp(:)~=statsSingleChannel.mask(:))
% figure
% imagesc(squeeze(statsStruct.mask(chech_ch,:,:)));
% figure
% imagesc(squeeze(statsSingleChannel.mask));

end


function  perform_TF_correction_nonparametric(tfr, pos)
% performs the cluster permutation test and uses plot2d for plotting in pos. of channels
% it is non-parametric which is suitable when the data is not guaranteed to
% be normal like accuracy ... it uses wilcoxon for the sample wise
% calculation of stats. between conditions ... no correction in space ..
% always assumes one channel

cond1_idx = mod(1:size(tfr.powspctrm,1),2);
cond1S =  tfr.powspctrm(logical(cond1_idx), :,:,:); % 'subj_chan_freq_time'
cond2S =  tfr.powspctrm(~logical(cond1_idx), :,:,:);

% cutting data if needed and setting chance level if any
chance = 0.5;  Tocut_temp=50;
% [chance, Tocut_temp] = lv_tune_params('if we have chance level type it','','if we need to cut a period from rows and cols (ms)',''); % chance for making limits of color bars dynamic, and cutting hint: half the smoothing window if erp features

Tocut = nearest(tfr.time, 0) + nearest(tfr.time, Tocut_temp/1000);
cond1S = cond1S(:,:,Tocut:end-Tocut,Tocut:end-Tocut); tfr.time = tfr.time(Tocut:end-Tocut); tfr.freq = tfr.freq(Tocut:end-Tocut);
cond2S = cond2S(:,:,Tocut:end-Tocut,Tocut:end-Tocut); tfr.powspctrm = tfr.powspctrm(:,:,Tocut:end-Tocut,Tocut:end-Tocut);


cfg = [];
cfg.test            = 'permutation';
cfg.correctm        = 'cluster';
cfg.n_permutations  = 1000;
cfg.clustercritval  = 1.96; % de el sample wise msh 3la el cluster el kbeer
% for alpha 0.05 = t-val is 1.96 ,,,, 1.645 for p=0.1
% for alpha 0.01 = t-val is 2.58
cfg.alpha = 0.05;
% Level 2 stats settings
cfg.design          = 'within';
cfg.statistic       = 'wilcoxon';
cfg.null            = 0 ; % the difference between conditions will be compared to 0

temp = tfr; temp.powspctrm=[];
% searching for cond1's positive effect only ...
for i=1:size(cond1S,2) % on every channel
    
    for subj = 1:size(cond1S,1)
        Diffresult{subj,1}.perf = squeeze(cond1S(subj,i,:,:)) - squeeze(cond2S(subj,i,:,:)); % .perf has freq_time
        Diffresult{subj,1}.metric = [];
    end
    
    statsStruct = mv_statistics(cfg, Diffresult);
    %statsStruct = lv_2d_clusterstats(cfg, Diffresult);
    % calculating the sample wise z-values
    x = squeeze(cond1S(:,i,:,:)); y = squeeze(cond2S(:,i,:,:)); zval=[];
    for k=1:size(cond1S,3) %freq
        for u=1:size(cond1S,4) %time
            if sum(isnan(squeeze(x(:,k,u))))==0
                [~,~,stats] = signrank(squeeze(x(:,k,u)),squeeze(y(:,k,u)), 'method','approximate');
                zval(k,u) = stats.zval; % freq_time
            else zval(k,u)=nan; end
        end
    end
    
    statsStruct.mask = double(statsStruct.mask);
    temp.powspctrm(i,:,:) = zval; % to make temp ch_freq_time
    mask = squeeze(statsStruct.mask);
    
    P_values = sort(statsStruct.p)
    % appending the channel label with the pVal of the significant clusters found
    Pid = unique(statsStruct.mask_with_cluster_numbers);
    Pid(Pid==0)=[];
    if ~isempty(Pid>0), temp.label = join([temp.label ', p= ' num2str(statsStruct.p(Pid))]); end
end

temp.cond1(size(tfr.powspctrm,2),:,:) = squeeze(mean(cond1S,1)); % ch_freq_time
temp.cond2(size(tfr.powspctrm,2),:,:) = squeeze(mean(cond2S,1));
temp.chance = chance;

plot2d(temp ,'z-stat, (cond1 positive clusters are shown)',pos, mask )


end



function perform_difference(data, pos)
warning('lv: data has one sbj we will get the difference between the two conditions');
cond1 = data.powspctrm(1, :,:,:); % sbj(cond)_ch_freq_time
cond2 = data.powspctrm(2, :,:,:);
% getting the difference
tfr_difference = data;
tfr_difference.powspctrm = squeeze( (cond1-cond2) ./ (cond1+cond2) ); % ch_freq_time

temp.powspctrm = tfr_difference.powspctrm; % was z-val but now difference becuase sbj lvl
temp.cond1 = cond1; % ch_freq_time
temp.cond2 = cond1;
temp.chance = 0; % assumed chance just to be able to plot the color coded erps !
plot2d(temp ,'Sbj lvl ratio: (cond1-cond2)./(cond1+cond2)',pos, [] );

end









function plot2d(data_struct ,colorlabel,pos, mask )
% takes data_struct with powspctrm ch_freq_time and the conditions as well and colorlabel and positions of channels
% then it loops on channels and plots in pos with unified colorbar scaling
% it plots the z-val or difference and the original conditions: cond1 and
% cond2
data = data_struct.powspctrm; % z-val or difference if one sbj
cond1 = data_struct.cond1; % getting conditions for plotting
cond2 = data_struct.cond2;
if isfield(data_struct,'chance')
    chance = data_struct.chance; end

if isnan(chance), x = [ max(cond1(:)) max(cond2(:)) min(cond1(:)) min(cond2(:)) ]; % if no chance then get the limits from data
    limits_conditions = [min(x) max(x)]; % non-symmetric
else
    x = max([max(cond1(:)) max(cond2(:))]); % symmetric around chance based on maximum
    limits_conditions = [chance-(x-chance) x];
end

x = max(abs(data(:))); % z-val has chance 0 .. so will be symmetric around that
limits = [-x x];


for i=1:size(data,1)
    if size(data,1)>1 % single figure if channels are compressed or if AUC, thus we only have one channel so we need one figure
        axes('Position',[pos(i,1) pos(i,2)  0.04 0.04]);
    end
    figure(100), % to plot on this specific figure I am giving it hypothetical ID
    b = imagesc(data_struct.time, data_struct.freq, squeeze(data(i,:,:)) ); set(gca,'YDir','normal'); %z-val plot, and unify the limits for all plots: caxis(limits)
    caxis(limits); title(data_struct.label(i)); set(b,'AlphaData',~isnan(squeeze(data(i,:,:)) ) ); % this AlphaData in case you want the cluster to be zval and anything else is white then make the mask ones and nans
    
    figure(200),
    imagesc(data_struct.time, data_struct.freq, squeeze(cond1(i,:,:)) ); set(gca,'YDir','normal'); %conditions plots
    caxis(limits_conditions);
    figure(300),
    imagesc(data_struct.time, data_struct.freq, squeeze(cond2(i,:,:)) ); set(gca,'YDir','normal');
    caxis(limits_conditions);
    
    figure(100),
    if ~isempty(mask) % if we have a mask then we plot its edges on top of the 2d figure of t-stat or z-stat
        BW = edge(mask);
        hold on,
        [row,col] = find(BW==1);
        plot(data_struct.time(col),data_struct.freq(row),'ZDataSource','','MarkerFaceColor',[1 1 1],...
            'MarkerEdgeColor',[1 1 1],...
            'Marker','square',...
            'LineStyle','none',...
            'MarkerSize',3,...
            'Color',[1 0 0]);
    end
end

if size(data,1)>1 % many channels and many figures
    figure(100), % the generalised color bar not the small ones near each small figure.. this one is the large one showing the unified color bar data.
    c = colorbar('peer',axes('Position',[0.88 0.1 0.1 0.8]));
    caxis(limits);
    c.Label.String = colorlabel; axis('off');  c.Label.Position(1) = -2;
    
    figure(200),
    c = colorbar('peer',axes('Position',[0.88 0.1 0.1 0.8]));
    caxis(limits_conditions);
    c.Label.String = colorlabel; axis('off');  c.Label.Position(1) = -2;
    figure(300),
    c = colorbar('peer',axes('Position',[0.88 0.1 0.1 0.8]));
    caxis(limits_conditions);
    c.Label.String = colorlabel; axis('off');  c.Label.Position(1) = -2;
else
    figure(100)
    h = colorbar; title(''); title(h, data_struct.label(i));
    xlabel('Wake time (sec.)'), ylabel('Sleep time (sec.)');
    figure(200)
    h = colorbar; title(h,'AUC'); title('Experimental night');
    xlabel('Wake time (sec.)'), ylabel('Sleep time (sec.)');
    figure(300)
    h = colorbar; title(h,'AUC'); title('Adaptation night');
    xlabel('Wake time (sec.)'), ylabel('Sleep time (sec.)');
end

end


 




% stat. helping functions
function statsStruct = lv_2d_clusterstats(cfg, result)
% 2d cluster permutation .. shuffle the sbj between cond and get the sample
% wise stats and sum the stats of contiguous points ... here we work on the
% difference between conditions so that if we want to shuffle a sbj
% between conditions we just flip the sign of the difference of the 'result'
% and compare to 0 and get the z-stats

% getting observed clusters
sample_alpha = cfg.clustercritval;
final_alpha = cfg.alpha;
permutations = cfg.n_permutations;

for i=1:length(result), cond_difference(i,:,:) = result{i,1}.perf; end
clear result;
tic
zval = zval_2d(cond_difference);
toc
L = bwlabel(zval>sample_alpha,8); % find the connected points in any direction '8-connectivity' and labels the connected clusters
fprintf(['\n Found: ' num2str(max(L(:))) ' clusters \n']);

for i=1:max(L(:)), ob_cluster_sum(i) = sum(sum(zval.*(L==i),1),2); end

decimal_val_lim = bi2de( ones(1,size(cond_difference,1)) ); % max decimal values if all sbj will be shuffled (all ones), one=flip the sign(shuffle)
shuffle_info = de2bi( randi(decimal_val_lim, [1,permutations]) );

clusters_counters = zeros(1,length(ob_cluster_sum));% counter of the times the shuffled data is higher than the observed val.
f = waitbar(0,'performing non-parametric cluster permutation');

% getting permutations
for i=1:permutations
    shuff_cluster_sum=[];
    cond_difference_shuffled = cond_difference;
    cond_difference_shuffled(shuffle_info(i,:)==1, :,:) = cond_difference(shuffle_info(i,:)==1, :,:) .* -1;
    zval_shuffled = zval_2d(cond_difference_shuffled);
    
    L_shuff = bwlabel(zval_shuffled>sample_alpha,8);
    for ii=1:max(L_shuff(:)), shuff_cluster_sum(ii) = sum(sum(zval_shuffled.*(L_shuff==ii),1),2); end
    clusters_counters( max(shuff_cluster_sum) > ob_cluster_sum ) = clusters_counters( max(shuff_cluster_sum) > ob_cluster_sum )+1;
    waitbar(i/permutations,f);
end
close(f);

statsStruct.p = clusters_counters ./ permutations;

statsStruct.mask = zeros(size(L));

idx = find(statsStruct.p<final_alpha);
for i=1:length(idx), statsStruct.mask(L==idx(i)) = 1; end % to be one despite the cluster

statsStruct.mask_with_cluster_numbers = statsStruct.mask .* L; % to show the actual cluster number

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






