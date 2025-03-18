function mo_pretty_errorbar(xax,cond1,cond2, stat,  varargin) %varargin for optional correction in space
% takes two conditions (sbj x time) and plots them across time with the standard error
% and performs point wise stats ( stat=0),
% ( stat=1) cluster stats with fieldtrip parametric t-test
% ( stat=2) cluster stats with matthias' non-parametric wilcoxon_test .. TODO because I need to see the inputs every time and also it works with his struct format

plot_individual=1; % plots individual curves with gray to be able to visually inspect outliers etc.
if plot_individual==1
    hold on, plot(xax,cond1,'Color',[0.878431379795074 0.878431379795074 0.878431379795074]);
    hold on, plot(xax,cond2,'Color',[0.541176497936249 0.831372559070587 0.686274528503418]);
end

if isempty(varargin)
    options=[];
    options.x_axis = xax;
    options.handle     = gcf;
    options.color_area = [128 193 219]./255;    % Blue theme
    options.color_line = [ 52 148 186]./255;
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'sem';
    plot_areaerrorbar(cond1,options )
    
    hold on;
    options=[];
    options.x_axis = xax;
    options.handle     = gcf;
    options.color_area = [243 169 114]./255;    % Orange theme
    options.color_line = [236 112  22]./255;
    options.alpha      = 0.5;
    options.line_width = 2;
    options.error      = 'sem';
    plot_areaerrorbar(cond2,options )
end



signi = nan(1,size(cond1,2));
if  stat==0 % significance is assumed both ways so either cond1 is higher or cond2 then it will be significant
    P_wilcoxon = nan(1,size(cond1,2));
    for i=1:size(cond1,2)
        [P_wilcoxon(1,i)] = signrank(cond1(:,i), cond2(:,i));
    end
    
    mask = zeros(1,size(cond1,2));
    mask(P_wilcoxon<0.05) = 1;
end

if  stat==1
    [ mask,~ ] = correct1D_Fieldtrip(xax,cond1,cond2,1:size(cond1,1), varargin);
    
end


if  stat==2
    [ mask,~ ] = perform_correction_nonparametric(xax,cond1,cond2,size(cond1,1));
    
end

if stat~=99
    hold on,
    axes_h = get(gcf,'CurrentAxes');
    signi(mask==1) = axes_h.YLim(2);
    signi(mask==0) = nan;
    % different color: [0.725490212440491 0.898039221763611 0.756862759590149]
    a = area(xax, signi, 'BaseValue',axes_h.YLim(1),'LineStyle','none', 'FaceColor',[0.235294118523598 0.831372559070587 0.0862745121121407]);
    a.FaceAlpha = 0.2;
    
    axis([axes_h.XLim(1) axes_h.XLim(2) axes_h.YLim(1) axes_h.YLim(2)]);
end


end





function [ mask,signi_time ] = correct1D_Fieldtrip(time_ax,cond1,cond2,sbj,  correct_in_space ) % the clustering window is the whole time window ... and two tailed by default
% format that matches: "load ERF_orig;" http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/#within-subjects-experiments
% formatting my data



for i=1:size(cond1,1), cond1S.trial{1,i} = cond1(i,:); end
for i=1:size(cond1,1), cond1S.time{1,i} = time_ax; end

cond1S.fsample = nearest(time_ax,time_ax(1)+1) - nearest(time_ax,time_ax(1));
warning(['sampling rate is set to: ' num2str(cond1S.fsample)]);
for i=1:size(cond2,1), cond2S.trial{1,i} = cond2(i,:); end
for i=1:size(cond2,1), cond2S.time{1,i} = time_ax; end
cond2S.fsample = cond1S.fsample;

for i=1:size(cond1S.trial,2)
    cond1Cells{1,i}.avg = cond1S.trial{1, i};  cond1Cells{1,i}.time = cond1S.time{1, 1};
    cond1Cells{1,i}.fsample = cond1S.fsample; if isempty(correct_in_space), cond1Cells{1,i}.label = {'measure-as-channel'}; end
    cond1Cells{1,i}.dimord = 'chan_time';
    
    cond2Cells{1,i}.avg = cond2S.trial{1, i};  cond2Cells{1,i}.time = cond2S.time{1, 1};
    cond2Cells{1,i}.fsample = cond2S.fsample; if isempty(correct_in_space), cond2Cells{1,i}.label = {'measure-as-channel'}; end
    cond2Cells{1,i}.dimord = 'chan_time';
end


%% actual test
cfg             = [];
cfg.method      = 'montecarlo'; % Monte Carlo approximation by creating the null hypothesis distribution
cfg.statistic   = 'depsamplesT';  % every UO is assigned to multiple experimental conditions in a particular order (within UO-design; dependent samples) [tha same subject in different conditions.. dependent]
cfg.clusteralpha = 0.05;        % alpha for thresholding the t statistic, a threshold ..

cfg.latency     =  [time_ax(1) time_ax(end)]; % clustering window

% Setting up the cluster test
cfg.correctm            = 'cluster'; % bonferoni, fdr etc.,
cfg.clusterstatistic    = 'maxsum';
cfg.alpha               = 0.05; % two-tailed =0.025 and change tail to 0 ... one tailed: I am testing en el exp 23la mn el adp msh en el mean bta3 el exp msh bysawe mean el adp ..
cfg.clustertail         = 1; % two tailed 0 ... 1 for one tailed
cfg.tail                = 1; % two tailed 0 ... 1 for one tailed
cfg.numrandomization    = 1000;

cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';

% Set up design matrix
% Row #1: subject nr (repeat for both conditions)
% Row #2: condition 1, then condition 2
design = zeros(2,2*numel(sbj));
design(1,:)= [1:numel(sbj), 1:numel(sbj)];
design(2,1:numel(sbj))        = 1;
design(2,numel(sbj)+1:2*numel(sbj)) = 2;

cfg.design      = design;
cfg.uvar        = 1;  % unit variable in 1st row
cfg.ivar        = 2;  % independent variable in 2nd row

if ~isempty(correct_in_space)&&correct_in_space==1
    load lv_layout lv_layout; cfg_temp = [];
    cfg_temp.method   = 'triangulation'; cfg_temp.senstype = 'EEG'; cfg_temp.layout = lv_layout;
    cfg.neighbours = ft_prepare_neighbours(cfg_temp);
    
    stat = ft_timelockstatistics(cfg, cond1Cells{:}, cond2Cells{:});
    cond1Cells.mask = stat.mask;
    
    cfg = [];
    cfg.layout = lv_layout;
    cfg.parameter = 'avg';
    cfg.maskparameter = 'mask';
    cfg.maskstyle = 'outline';
    ft_multiplotER(cfg, cond1Cells);
end


stat = ft_timelockstatistics(cfg, cond1Cells{:}, cond2Cells{:});


if (isfield(stat,'posclusterslabelmat')==1) && length(stat.posclusters)>0
    %hold on,
    %plotting positive cluster(s)
    %plot(stat.time, (stat.posclusterslabelmat).*(stat.mask));
    fprintf('positive clusters Ps:')
    Ps = squeeze( struct2cell(stat.posclusters ) )';
    (cell2mat(Ps(:,1)))
    
    maxStat = max(cell2mat(Ps(:,2)))
end
if (isfield(stat,'negclusterslabelmat')==1) && length(stat.negclusters)>0
    %     hold on,
    %plotting negative cluster(s)
    %     plot(stat.time, (stat.negclusterslabelmat).*(stat.mask));
    fprintf('min negative clusters Ps:')
    Ps = squeeze( struct2cell(stat.negclusters ) )';
    min(cell2mat(Ps(:,1)))
end
mask = stat.mask;

signi_time = stat.time(mask);
end

%% based on Matthias MVPA non-parametric correction ..
function [ mask,signi_time ] =  perform_correction_nonparametric(time_ax,cond1,cond2,sbj,  correct_in_space )
% performs the cluster permutation test for curves
% it is non-parametric which is suitable when the data is not guaranteed to
% be normal like accuracy ... it uses wilcoxon for the sample wise
% calculation of stats. between conditions ... no correction in space ..
% always assumes one channel

cond1S =  cond1; % 'subj_time'
cond2S =  cond2;

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
cfg.null            = 0; % the difference between conditions will be compared to 0

% searching for cond1's positive effect only ...
for subj = 1:size(cond1S,1)
    Diffresult{subj,1}.perf = (squeeze(cond1S(subj,:)) - squeeze(cond2S(subj,:)))'; % .perf has freq_time
    Diffresult{subj,1}.metric = [];
end

warning('using lv_cluster_permutation, to use MVPA uncomment MVPA stats part in lv_plot_topo.m')
[stat,~] = lv_mdim_clusterstats(cfg, Diffresult, 0);

% MVPA stats 
% stat = mv_statistics(cfg, Diffresult);



fprintf(['\n This is a positive clusters only test we found: P = ' num2str(stat.p) '\n']);

mask = stat.mask;

signi_time = time_ax .* double(stat.mask);

end






