function [ erps_struct ] = lv_erp(varargin) %data, do_stats, do_plot
% calculates the ERPs for single participant or group level and performs stats. on the
% result by calling lv_plot_topo
% if group level then data.trial should contain erps of cond1 then cond2 on
% the first dim so the odd idx is cond1 and cond2 is even

fprintf(['\n Performing ERP analysis \n']);

if nargin < 3  do_plot = 0; end
if nargin < 2  do_plot = 0; do_stats = 0; end
if varargin{2}==1  do_plot = 1; end

data = varargin{1};
do_stats = varargin{2};
do_plot = varargin{3};
load lv_layout lv_layout;

% not 3d make it 3d
if length(size(data.trial))<3
    cfg             = [];
    cfg.keeptrials  = 'yes';
    data = ft_timelockanalysis(cfg, data);
end

% interpolating channels locations to fit in many subplots on the screen
pos = [interp1([min(lv_layout.pos(:,1)) max(lv_layout.pos(:,1))],[0.1 0.9],lv_layout.pos(:,1)) ...
    interp1([min(lv_layout.pos(:,2)) max(lv_layout.pos(:,2))],[0.1 0.9],lv_layout.pos(:,2))];

% if there is trialinfo then we have one sbj, level1
if isfield(data,'trialinfo')
    % erps calculation
    conds = unique(data.trialinfo(:,1)); if length(conds)>2, error('data has more than two conditions!'); end
    
    baseline = [ ]; % ex: [-0.5 0]
    fprintf(['baseline is set to be: ' num2str(baseline) '\n']);
    if ~isempty(baseline)
        % single trial baseline
        baselines_mx = mean(data.trial(:,:,nearest(data.time,baseline(1)):nearest(data.time,baseline(2))), 3); % trl_ch
        baselines_mx = repmat(baselines_mx,[1 1 length(data.time)]);
        data.trial = data.trial - baselines_mx;
    end
    erps(1,:,:) = mean(data.trial(data.trialinfo(:,1)==conds(1),:,:),1);
    erps(2,:,:) = mean(data.trial(data.trialinfo(:,1)==conds(2),:,:),1);
    
    do_stats=0; % for sbj level only visualise don't do stats because for sbj lvl we can't do it in the same
    % way as the group because the no. trials for conditions is
    % different so we can't calculate the stats.. and trials are independent
    if do_plot==1
        lv_plot_topo(data, conds, do_stats, pos);
    end
    
    erps_struct = data;
    erps_struct = rmfield(erps_struct,'trialinfo');
    erps_struct.trial = erps;
    return;
else % group level
    data.lv_layout = lv_layout;
    lv_plot_topo(data, [], do_stats, pos);
    
end



end









