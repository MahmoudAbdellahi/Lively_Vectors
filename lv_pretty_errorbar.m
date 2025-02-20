function all_stats = lv_pretty_errorbar(varargin)
% performs stats between vectors of conditions .. either vector vs vector
% when you give two vectors only .. or across time with correction in time
% and also optional correction in space and also parametric or
% non-parametric or point-wise stats
 
% because when across time the first arg. is timeax which equals the second
% dim. of the cond1(sbj_time)
if numel(varargin{1})==size(varargin{2},2)
    xax=varargin{1}; cond1=varargin{2}; cond2=varargin{3};
    stat=varargin{4}; correct_in_space=[]; if length(varargin)==5, correct_in_space=varargin{5}; end
%     figure,
    all_stats = pretty_errorbar_core(xax,cond1,cond2, stat,  correct_in_space);
else
    t1=varargin{1}; t2=varargin{2};
    
    [ all_stats ] =  vec2stat(t1,t2)
end


end




%% core functions
%% two vectors
function [ all_stats ] =  vec2stat(t1,t2)
% takes two vectors and calculates the tstat and wilcoxon
% also plots them... if t2 is chance then we plot t1
figure,
all_stats.P_wilcoxon=nan;

if length(t1)==length(t2) % paired samples with equal lengths
    [~,all_stats.P_t_test,all_stats.confidenceInterval,all_stats.Tstats] = ttest(t1, t2);
    [all_stats.P_wilcoxon,~,all_stats.Wilcoxon_stats] = signrank(t1, t2, 'method' ,'approximate');
    all_stats.Tstats = all_stats.Tstats.tstat;
else
    [all_stats.P_wilcoxon,~,all_stats.Wilcoxon_stats] = ranksum(t1, t2, 'method' ,'approximate');
    warning('comparing to one value and using ranksum not signrank ....')
end
all_stats.mean1 = mean(t1);
all_stats.mean2 = mean(t2);
all_stats.wilcoxonZval = all_stats.Wilcoxon_stats.zval;

r = linspace(1,1.1,length(t1));
r2 = linspace(2,2.1,length(t2));


errorbar(0.95, mean(t1,1) , std(t1,[],1)./sqrt(size(t1,1)),'Marker','square',...
    'LineWidth',2,'color', [0 0.4470 0.7410]),
hold on, scatter(r, t1, 50, 'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532], 'LineWidth',0.75);
if length(t2)>1
    hold on,
    errorbar(1.95, mean(t2,1) , std(t2,[],1)./sqrt(size(t2,1)),'Marker','square',...
        'LineWidth',2, 'color',[0.8500 0.3250 0.0980]),
    hold on, scatter(r2, t2,50, 'MarkerEdgeColor',[0.850980401039124 0.325490206480026 0.0980392172932625], 'LineWidth',0.75);
    
    [xlb,ylb,connect_pts,chance] = lv_tune_params('xlabel','Cond1','ylabel','Cond2','connect conditions pts?','0','chance level','');
    
    h = gca;
    set(h, 'FontSize',18,'TickLabelInterpreter','none','XTick',[1 2],...
        'XTickLabel',{string(xlb),string(ylb)}, 'box','off');
    
    
    xlim([0.5  2.5])
    
    if connect_pts==1
        for i=1:length(t1)
           hold on, 
           plot([r(i) r2(i)],[t1(i) t2(i)],'color',[0.80,0.80,0.80]); 
        end
    end
    if isnumeric(chance)
       hold on, plot([0.5 2.5],[chance chance],'LineStyle','--','color',[0.80,0.80,0.80]); 
    end
end
if length(t2)==1
    hold on, plot([0.5 1.5],[t2 t2]);
    
    [xlb] = lv_tune_params('xlabel','Cond1');
    
    h = gca;
    set(h, 'FontSize',18,'TickLabelInterpreter','none','XTick',[1],...
        'XTickLabel',{string(xlb)}, 'box','off');
     
    xlim([0.5  1.5])
end
 
title(['Pwilcoxon=' num2str(all_stats.P_wilcoxon)],'FontSize',11);
end

%% across time 
function stats=pretty_errorbar_core(xax,cond1,cond2, stat,  varargin) %varargin for optional correction in space
% takes two conditions (sbj x time) and plots them across time with the standard error
% and performs point wise stats ( stat=0),
% ( stat=1) cluster stats with fieldtrip parametric t-test
% ( stat=2) cluster stats with matthias' non-parametric wilcoxon_test+
p_vals = nan; mask=[];
plot_individual=0; % plots individual curves with gray to be able to visually inspect outliers etc.
if plot_individual==1
    hold on, plot(xax,cond1,'Color',[0.878431379795074 0.878431379795074 0.878431379795074]);
    hold on, plot(xax,cond2,'Color',[0.541176497936249 0.831372559070587 0.686274528503418]);
end

if  stat==0 % significance is assumed both ways so either cond1 is higher or cond2 then it will be significant
    P_wilcoxon = nan(1,size(cond1,2));  %r = (xax(2)-xax(1))/1.5; % place of errorbar .. use it for log spacing 
    for i=1:size(cond1,2)
        hold on, %scatter(xax(i), cond1(:,i), 50, 'MarkerEdgeColor',[0.62,0.77,0.88], 'LineWidth',0.25);
        [P_wilcoxon(1,i)] = signrank(cond1(:,i), cond2(:,i)); 
        errorbar(xax(i), mean(cond1(:,i),1) , std(cond1(:,i),[],1)./sqrt(size(cond1(:,i),1)),'Marker','square',...
        'LineWidth',2,'color', [0,0.45,0.74]);
        errorbar(xax(i), mean(cond2(:,i),1) , std(cond2(:,i),[],1)./sqrt(size(cond2(:,i),1)),'Marker','square',...
        'LineWidth',2,'color', [1.00,0.41,0.16]);
    end 
    for i=1:length(P_wilcoxon)
        r = ylim;
        if P_wilcoxon(1,i)<0.05, plot(xax(i),r(2),'*','color',[0.15,0.15,0.15]); end
        P_wilcoxon(1,i)
    end
    h = gca;
    set(h, 'FontSize',10,'box','off');
    xticks(xax), xax(P_wilcoxon<0.05)
end

varargin=cell2mat(varargin);
if isempty(varargin) && stat~=0
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


if  stat==1
    [ mask,~,p_vals,maxStat ] = perform_parametric_correction(xax,cond1,cond2,1:size(cond1,1), varargin);
     stats.maxStat=maxStat; 
end


if  stat==2
    [ mask,~,p_vals ] = perform_nonparametric_correction(xax,cond1,cond2,size(cond1,1));
   
end

% for lv, 11 for parametric and 22 for nonparametric
if stat==11
    parametric = 1; [ mask,~,p_vals ] = perform_correction_with_lv(parametric, xax,cond1,cond2,1:size(cond1,1), varargin); 
end
if stat==22
    parametric = 0; [ mask,~,p_vals ] = perform_correction_with_lv(parametric, xax,cond1,cond2,1:size(cond1,1), varargin); 
end



if stat~=99 && stat~=0
    hold on,
    axes_h = get(gcf,'CurrentAxes');
    signi(mask==1) = axes_h.YLim(2);
    signi(mask==0) = nan;
    % different color: [0.725490212440491 0.898039221763611 0.756862759590149]
    a = area(xax, signi, 'BaseValue',axes_h.YLim(1),'LineStyle','none', 'FaceColor',[0.235294118523598 0.831372559070587 0.0862745121121407]);
    a.FaceAlpha = 0.2;
    
    axis([axes_h.XLim(1) axes_h.XLim(2) axes_h.YLim(1) axes_h.YLim(2)]);
end

title(['P= ' num2str(p_vals(:)')],'FontSize',11);

stats.mask=mask;
stats.p_vals=p_vals;

end





function [ mask,signi_time,p_vals,maxStat ] = perform_parametric_correction(time_ax,cond1,cond2,sbj,  correct_in_space ) % the clustering window is the whole time window
% format that matches: "load ERF_orig;" http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/#within-subjects-experiments
% formatting my data
correct_in_space = cell2mat(correct_in_space);


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
cfg.statistic   = 'depsamplesT';  % every UO is assigned to multiple experimental conditions in a particular order (within UO-design; dependent samples) [tha same participant in different conditions.. dependent]
cfg.clusteralpha = 0.05; % alpha for thresholding the t statistic

cfg.latency     =  [time_ax(1) time_ax(end)]; % clustering window

% Setting up the cluster test  
% [cfg.correctm] = lv_tune_params('do you want cluster(0) or TFCE (1) classification?','0');
cfg.correctm=0;
if cfg.correctm==0, cfg.correctm = 'cluster'; % bonferoni, fdr etc.,
    cfg.clusterstatistic    = 'maxsum';
    cfg.alpha               = 0.05; 
    cfg.tail                = 1; 
    cfg.clustertail         = 1; 
    cfg.numrandomization    = 10000;
else
    cfg.correctm = 'tfce'; cfg.tfce_H = 2; cfg.tfce_E = 0.5; % default setting
    cfg.numrandomization = 'all';   % there are X participants, so 2^X possible raondomizations
end 

cfg.avgoverchan = 'no';
cfg.avgovertime = 'no';

design = zeros(2,2*numel(sbj));
design(1,:)= [1:numel(sbj), 1:numel(sbj)];
design(2,1:numel(sbj))        = 1;
design(2,numel(sbj)+1:2*numel(sbj)) = 2;

cfg.design      = design;
cfg.uvar        = 1;  % unit variable
cfg.ivar        = 2;  % independent variable


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
  
p_vals=nan; maxStat=nan;
if isfield(stat,'posclusters')
    maxStat = cell2mat(squeeze( struct2cell(stat.posclusters ) ));
    maxStat = max(maxStat,[],2); maxStat = maxStat(2);
end
if (isfield(stat,'posclusterslabelmat')==1) && length(stat.posclusters)>0
    %hold on,
    %plotting positive cluster(s)
    %plot(stat.time, (stat.posclusterslabelmat).*(stat.mask));
    fprintf('positive clusters Ps:')
    Ps = squeeze( struct2cell(stat.posclusters ) )';
    p_vals = (cell2mat(Ps(:,1)))
    
    %maxStat = max(cell2mat(Ps(:,2)))
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
function [ mask,signi_time,p_vals ] =  perform_nonparametric_correction(time_ax,cond1,cond2,sbj,  correct_in_space )
% performs the cluster permutation test for curves
% it is non-parametric which is suitable when the data is not guaranteed to
% be normal like accuracy ... it uses wilcoxon for the sample wise
% calculation of stats. between conditions ... no correction in space ..

cond1S =  cond1; % 'ppnt_time'
cond2S =  cond2;

cfg = [];
cfg.test            = 'permutation';
cfg.correctm        = 'cluster';
cfg.n_permutations  = 10000;
cfg.clustercritval  = 1.96; % sample alpha
cfg.alpha = 0.05;
% Level 2 stats settings
cfg.design          = 'within';
cfg.statistic       = 'wilcoxon';
cfg.null            = 0; % the difference between conditions vs. 0
 
% searching for cond1's positive effect
for subj = 1:size(cond1S,1)
    Diffresult{subj,1}.perf = (squeeze(cond1S(subj,:)) - squeeze(cond2S(subj,:)))';
    Diffresult{subj,1}.metric = [];
end

stat = mv_statistics(cfg, Diffresult);

fprintf(['\n Checked for positive clusters and found: P = ' num2str(stat.p) '\n']);

mask = stat.mask;

signi_time = time_ax .* double(stat.mask);
p_vals=nan;
if isfield(stat,'p'), p_vals = stat.p; end

end

%% LV correction
function [ mask,signi_time,p_vals ] =  perform_correction_with_lv(parametric, time_ax,cond1,cond2,sbj,  correct_in_space )
% performs the cluster permutation test for curves
% with lv, parametric or nonparametric

cond1S =  cond1; % 'ppnt_time'
cond2S =  cond2;

cfg = [];
cfg.test            = 'permutation';
cfg.correctm        = 'cluster';
cfg.n_permutations  = 100000;
cfg.clustercritval  = 1.96; % sample alpha 
cfg.alpha = 0.05;
% Level 2 stats settings
cfg.design          = 'within';
cfg.statistic       = 'wilcoxon';
cfg.null            = 0; % the difference between conditions vs. 0
 
% searching for cond1's positive effect
for subj = 1:size(cond1S,1)
    Diffresult{subj,1}.perf = (squeeze(cond1S(subj,:)) - squeeze(cond2S(subj,:)))';
    Diffresult{subj,1}.metric = [];
end
 
stat = lv_mdim_clusterstats(cfg, Diffresult, parametric);

fprintf(['\n Checked for positive clusters and found: P = ' num2str(stat.p) '\n']);

mask = stat.mask;

signi_time = time_ax .* double(stat.mask);
p_vals=nan;
if isfield(stat,'p'), p_vals = stat.p; end

end






