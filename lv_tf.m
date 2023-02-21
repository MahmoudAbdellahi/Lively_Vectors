function [ TF_struct ] = lv_tf(varargin) %data, do_stats, do_plot
% calculates the TF analysis for single sbj or group level and does stats on the
% result and plots the clusters

% if you have AUC or any other 2d data .. input:cond_ch_yaxis_xaxis (for txt classification:cond_ch_trn_tst)
% % EXAMPLE for txt classsification:
% % dat=[];
% % for i=1:size(result_exp,1), cond1(i,:,:)=result_exp{i,1}.perf; cond2(i,:,:)=result_adp{i,1}.perf; end
% % id = mod(1: size(cond1,1)+size(cond2,1), 2);
% % dat.trial(id==0,1,:,:) = cond2; % chance: (cond2.*0)+0.5;
% % dat.trial(id==1,1,:,:) = cond1;
% % dat.time=times.timeTst; dat.freq=times.timeTrn; dat.label={'z-stat'};
% % [ TF_struct ] = lv_tf(dat,1,1); %data, do_stats, do_plot

% gives rpt_ch_freq_time ... if sbj lvl then ch_freq_time
% and if group it gives rpt_ch_freq_time .. but rpt has cond1 then cond2
% for a particular sbj



fprintf(['\n Performing time frequency analysis \n']);

if nargin < 3  do_plot = 0; end
if nargin < 2  do_plot = 0; do_stats = 0; end
if varargin{2}==1  do_plot = 1; end

data = varargin{1};
do_stats = varargin{2};
do_plot = varargin{3};

if isfield(data,'lv_layout')
    lv_layout = data.lv_layout;
else
    load lv_layout lv_layout;
end
% not 3d make it 3d
if length(size(data.trial))<3
    cfg             = [];
    cfg.keeptrials  = 'yes';
    data = ft_timelockanalysis(cfg, data);
end

% interpolating channels locations to fit in many subplots on the screen
cfg = [];
cfg.layout = 'easycapM1.mat'; % will get just the max and minimum coordinates to plot in logical center positions
reference_lay = ft_prepare_layout(cfg);

pos = [interp1([min(reference_lay.pos(:,1)) max(reference_lay.pos(:,1))],[0.1 0.9],lv_layout.pos(:,1)) ...
    interp1([min(reference_lay.pos(:,2)) max(reference_lay.pos(:,2))],[0.1 0.9],lv_layout.pos(:,2))];


% if there is trialinfo then we are in one sbj, level1
if isfield(data,'trialinfo')
    % TF calculation
    conds = unique(data.trialinfo(:,1)); if length(conds)>2, error('data has more than two conditions!!'); end
    
    if isfield(data,'method')
        if strcmp(data.method,'itpc'), [TF(1,:,:,:), TFdat] = do_tf(data, [data.time(1) data.time(end)], [], [0.5 25], data.method); end % dat, window, baseline, frequencies, method
    else
        cfg= []; cfg.trials = find(data.trialinfo(:,1)==conds(1));
        [TF(1,:,:,:), TFdat] = do_tf(ft_selectdata(cfg, data), [data.time(1) data.time(end)], [-0.3 -0.1], [1 30], []); % dat, window, baseline, frequencies, method
        
        cfg= []; cfg.trials = find(data.trialinfo(:,1)==conds(2));
        [TF(2,:,:,:), ~] = do_tf(ft_selectdata(cfg, data), [data.time(1) data.time(end)], [-0.3 -0.1], [1 30], []);
    end
    
    do_stats=0; % for sbj level only visualise don't do stats
    
    TFdat.trial = TF; % just for .trial check in lv_plot_topo
    TFdat.powspctrm = TF;
    if do_plot==1
        lv_plot_topo(TFdat, [], do_stats, pos); % data, conds, do_stats, pos
    end
    
    TF_struct = TFdat;
    return;
else % group level
    data.powspctrm = data.trial;
    data.lv_layout = lv_layout;
    lv_plot_topo(data, [], do_stats, pos);
    TF_struct=[];
end



end




function [fullpow , TFdat]= do_tf(dat, window, baseline, frequencies, method) % takes 3d in .trial (trls_ch_time) and returns (ch_freq_time)

cfg              = [];
cfg.output       = 'pow';
if strcmp(method,'itpc')
    cfg.output       = 'fourier'; % to get the complex signal and extract the phase, mtmconvol is still performing a wavelet but fourier in the o/p
    % means you get the complex analytical signal
end
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = linspace(frequencies(1),frequencies(end),2*(1+frequencies(end)-frequencies(1)));% frequencies(1):0.5:frequencies(end) default: 2 to 30 HZ
% foi is arbitrary if big then it will repeat values but still correct
cfg.t_ftimwin    = 5./cfg.foi;  %default: 5 cycles per time window as a minimum to describe the frequency well
cfg.toi          = dat.time; % .time for max resolution .. window(1):0.1:window(2) the jumps just for visual smoothing
cfg.pad          ='nextpow2'; % faster estimation
TFdat = ft_freqanalysis(cfg, dat);
TFdat.fourierspctrm = single(TFdat.fourierspctrm);
if strcmp(method,'itpc')
    F = TFdat.fourierspctrm; TFdat=[];   % copy the Fourier spectrum
    N = size(F,1);                       % number of trials
    TFdat.powspctrm = F./abs(F);         % divide by magnitudes... to make all the magnitudes ones
    TFdat.powspctrm = sum(TFdat.powspctrm,1);   % sum angles .. because they are already to the exp power so just sum .. they already in euler's form
    TFdat.powspctrm = abs(TFdat.powspctrm)/N;   % take the absolute value .. because the magnitude reflects the coherence
    TFdat.powspctrm = squeeze(TFdat.powspctrm); % remove the first singleton dimension .. that was trials
end

if ~isempty(baseline)
    cfg              = [];
    cfg.baseline     = [baseline(1) baseline(2)];
    cfg.baselinetype = 'relchange';
    [TFdat] = ft_freqbaseline(cfg, TFdat); % ch x freq x time
end

fullpow = TFdat.powspctrm;
% ft_singleplotTFR(cfg, TFdat);


end







