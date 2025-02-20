%% Targeting targeted memory reactivation: characteristics of cued reactivation in sleep
% Mahmoud E. A. Abdellahi, Anne C. M. Koopman, Matthias S. Treder & Penelope A. Lewis

% This script is designed to compute temporal measures
% by Miguel Navarrete
% CUBRIC
% 2017
function varargout = lv_SO_Spindle_FeatureExtraction(correctness_trls, spindle)
%% Formatting data for the extraction of SO and spindle features
for i=1:length(correctness_trls)
    fs=200;
    new_header{i,1}.fsample = fs;
    st = ceil(size(correctness_trls{i, 1}.trial,3)/2);
    new_header{i,1}.label = correctness_trls{i, 1}.label;
    new_header{i,1}.sampleinfo = [st   size(correctness_trls{i, 1}.trial,1)*size(correctness_trls{i, 1}.trial,3) + st];
    new_header{i,1}.scoring = {  int32( ones( 1,  size(correctness_trls{i, 1}.trial,1)*size(correctness_trls{i, 1}.trial,3) ) )  };
    the1_2_stim = repmat([1;2], ceil( length(correctness_trls{i, 1}.sampleinfo(:,1)))  ,1);
    the1_2_stim = the1_2_stim(1:length(correctness_trls{i, 1}.sampleinfo(:,1)) *2);

    idx = st: size(correctness_trls{i, 1}.trial,3): new_header{i,1}.sampleinfo(2)+ st;
    dummy = idx ; % generating dummy idx for second click

    idx = idx(1:sum(the1_2_stim==1));
    dummy = dummy(1:sum(the1_2_stim==2));
    V = [];
    V(the1_2_stim==1) = idx;
    V(the1_2_stim==2) = dummy + 301; % to cut the samples according to the
    % trial length (1.5sec. in sleep and sampling rate is 200samples/sec., so 301samples)
    new_header{i,1}.stim  = [V'   the1_2_stim  ];
    new_header{i,1}.time{1, 1} = 0 : 1/fs : (size(correctness_trls{i, 1}.trial,1)*size(correctness_trls{i, 1}.trial,3)  )/200;
    new_header{i,1}.time{1, 1} = new_header{i,1}.time{1, 1}(1:size(correctness_trls{i, 1}.trial,1)*size(correctness_trls{i, 1}.trial,3));

    temp = permute(correctness_trls{i, 1}.trial,[2 3 1]);
    temp = reshape(temp,length(correctness_trls{i, 1}.label),[] );
    new_header{i,1}.trial{1, 1} = temp;
end 
%% Setting parameters
sbj = 1:length(correctness_trls); % codes of participants
channel = find(ismember(lower(new_header{1, 1}.label),lower({'Fz'})));
nm_fsample = 200;
motor = find(ismember(lower(new_header{1, 1}.label),lower({'C6','CP4','C5','CP3'}))); % motor 

vt_timeLim  = [-3.04, 3.04]; % long trial duration inside sleep_raw
nm_SMOOTHWINDOW     = 0.2; % 200ms window for rms power calculation
nm_TROUGHTHRESHOLD	= 0;  
nm_PEAKTHRESHOLD	= 0; % considering all SOs, for now and filter them later inside SO_Spindle classification
nm_STIMDISTANCE     = 2; % pre-cue duration for feature extraction
if spindle ==1, nm_STIMDISTANCE = 1.5; end % pre-cue duration for feature extraction

% detection thresholds
vt_timeSpindles     = [0.3,3];
vt_timeTroughs      = [0.166667,1.666667];
vt_timeTroughs      = vt_timeTroughs + [-1,+1].*diff(vt_timeTroughs)*0.05; % +/- 5% time range
vt_stageSpindles    = [3,4];
vt_stageSO          = [2,3,4];
vt_freqSpindleSlow	= [9,11];
vt_freqSpindleFast	= [11,17];
nm_SPINDLENUMOSC	= 5;


vt_timeEpoch        = vt_timeLim(1):1/nm_fsample:vt_timeLim(2);
ob_filtEEGLo= fn_designDCfirfilter(nm_fsample);
ob_filtEEGHi= fn_designLowPassEEGfirfilter(nm_fsample);
ob_filtFS	= fn_designSpindlefirfilter(nm_fsample,'custom',vt_freqSpindleFast);
ob_filtSS	= fn_designSpindlefirfilter(nm_fsample,'custom',vt_freqSpindleSlow);
ob_filtSO	= fn_designSOfirfilter(nm_fsample);

vt_hilFreq  = (0:0.01:nm_fsample/2);

vt_idxAmp1   = vt_hilFreq < 3;
vt_idxAmp2   = vt_hilFreq >= 3 & vt_hilFreq <= 4;

vt_idxAmp3   = vt_hilFreq > 4;

vt_idxAmp1   = ones(sum(vt_idxAmp1),1);
vt_idxAmp2   = (1:sum(vt_idxAmp2)).*(-1/sum(vt_idxAmp2))+1;
vt_idxAmp3   = zeros(sum(vt_idxAmp3),1);

vt_hilAmpl  = vertcat(vt_idxAmp1,vt_idxAmp2(:),vt_idxAmp3);
vt_hilFreq  = vt_hilFreq * 2 /nm_fsample;

if mod(numel(vt_hilFreq),2)
    vt_hilFreq  = vt_hilFreq(1:end-1);
    vt_hilAmpl  = vt_hilAmpl(1:end-1);
end

nm_hilOrder	= nm_fsample/2;
vt_filtImag	= firls(nm_hilOrder,vt_hilFreq,vt_hilAmpl,'hilbert');
vt_filtReal	= firls(nm_hilOrder,vt_hilFreq,ones(size(vt_hilAmpl)));
nm_hilDelay	= nm_hilOrder/2;

nm_smoothSamples    = round(nm_SMOOTHWINDOW * nm_fsample);


vt_dataLabels	= {'markerTime','trialsEEG','trialsSO','trialsFS','trialsSS'...
    'trialsFSrms','trialsSSrms','trialsSOphH','trialsSSphF'...
    'features','detectionLogs'};

%% Feature exraction
for ff = 1 : numel(sbj)
    fprintf('\n\n\n ---*** Processing participant %d ***--- \n\n', sbj(ff))

    st_Header = new_header{ff,1};
    vt_channels     = st_Header.label;
    st_Header.scoring{1, 1} = st_Header.scoring{1, 1}.*3; % all data from SWS

    %% - Set signals
    % Selcet & Obtain raw channels
    if numel(st_Header.label) == 1
        st_Header.label = st_Header.label{1};
    end
    mx_rawEEG   = nan(size(st_Header.trial{1},2),numel(vt_channels));
    for cc = 1:numel(vt_channels)
        vt_idxCh	=  ismember(st_Header.label,vt_channels{cc});

        if any(vt_idxCh)
            vt_tempChannel  = st_Header.trial{1}(vt_idxCh,:);
            mx_rawEEG(:,cc)	= vt_tempChannel(:);
        else
            warning('Channel %s not present ',...
                vt_channels{cc})
            mx_rawEEG(:,cc)	= 0;
        end
        clear vt_tempChannel
    end

    % Obtain sampling frequency
    nm_rawSampling	= st_Header.fsample;

    % Obtain markers
    mx_markers  = st_Header.stim;

    vt_idxStim                  = mx_markers(:,2) == 17;
    mx_markers(vt_idxStim,2)	= 1;
    vt_idxStim                  = mx_markers(:,2) == 18;
    mx_markers(vt_idxStim,2)	= 2;


    vt_idxStim  = mx_markers(:,2) == 1 | mx_markers(:,2) == 2;
    mx_markers  = mx_markers(vt_idxStim,:);

    % Prepare markers by trials
    mx_markers	= fn_preparemarkers(mx_markers,[1 2]);
    mx_markTime	= mx_markers./nm_rawSampling;

    % Obtain hypnogram
    vt_hypnogram	= st_Header.scoring{1};
    % Clear memory
    clear st_Header
    %% - Resample Data
    fprintf('\n Resampling data to %iHz: ',nm_fsample)
    mx_rawEEG	= resample(mx_rawEEG,nm_fsample,nm_rawSampling);

    vt_rawTime  = (0:numel(vt_hypnogram)-1)/nm_rawSampling;
    vt_eegTime  = (0:size(mx_rawEEG,1)-1)/nm_fsample;

    vt_hypnogram	= single(vt_hypnogram);
    vt_hypnogram	= interp1(vt_rawTime,vt_hypnogram,vt_eegTime,'nearest');
    vt_hypnogram	= int32(vt_hypnogram);
    %% - Process channels
    mx_dataClick_1	= cell(numel(vt_channels),numel(vt_dataLabels));
    mx_dataClick_2  = cell(numel(vt_channels),numel(vt_dataLabels));

    mx_spindleDet   = cell(numel(vt_channels),2);
    vt_soWaves      = cell(numel(vt_channels),1);
    vt_soAmpli      = cell(numel(vt_channels),1);
    vt_soThres      = cell(numel(vt_channels),1);
    vt_pWelch       = cell(numel(vt_channels),1);

    if spindle ==1, loopidx = motor; end
    if spindle ==0, loopidx = channel; end

    for cci = 1:length(loopidx ) % looping on channels
        cc = loopidx(cci);
        if sum(abs(mx_rawEEG(:,cc))) == 0
            continue
        end

        fprintf('\n	** Filtering in EEG band | ch-%s: ',...
            vt_channels{cc})

        vt_signalEEG	= fn_filterOffline(mx_rawEEG(:,cc),ob_filtEEGLo);
        vt_signalEEG	= fn_filterOffline(vt_signalEEG,ob_filtEEGHi);

        [st_psd.power,st_psd.freq]	= pwelch(vt_signalEEG,...
            10*nm_fsample,2*nm_fsample,...
            [],nm_fsample);

        vt_pWelch{cc}	= st_psd;

        % Filter in the fast spindle band
        vt_signalFS     = fn_filterOffline(mx_rawEEG(:,cc),ob_filtFS);
        vt_rmsFS        = fn_rmstimeseries(vt_signalFS,nm_smoothSamples);

        % Filter in the slow spindle band
        vt_signalSS     = fn_filterOffline(mx_rawEEG(:,cc),ob_filtSS);
        vt_rmsSS        = fn_rmstimeseries(vt_signalSS,nm_smoothSamples);

        % Filter in SO frequency band
        vt_signalSO     = fn_filterOffline(mx_rawEEG(:,cc),ob_filtSO);

        % Apply Hilbert filter
        vt_Imag	= filter(vt_filtImag,1,vt_signalSO);
        vt_Real	= filter(vt_filtReal,1,vt_signalSO);

        vt_hilb = vt_Real + 1i*vt_Imag;

        vt_hilb(1:end-nm_hilDelay)	= vt_hilb(nm_hilDelay+1:end);
        vt_hilb(end-nm_hilDelay:end)= nan;

        vt_hilbFilter	= vt_hilb;
        % Obtain SO phase
        vt_phaseSO	= angle(hilbert(vt_signalSO));

        %% Detect SO events
        fprintf('\n   ** Detect SO events: ')
        st_Cnf              = struct;
        st_Cnf.freqband     = [];
        st_Cnf.fsampling	= nm_fsample;
        st_Cnf.threshold	= nm_TROUGHTHRESHOLD;
        st_Cnf.timebounds   = vt_timeTroughs;
        st_Cnf.hypnogram	= vt_hypnogram;
        st_Cnf.method       = 'threshold';
        st_Cnf.stage        = vt_stageSO;
        st_Cnf.minthresh	= [];
        st_Cnf.toFilter     = [];

        [vt_SOwaves,nm_thr]	= fn_detectsleepSO(vt_signalSO,st_Cnf);
        vt_soWaves{cc}      = single(vt_SOwaves);
        vt_soAmpli{cc}      = single(vt_signalSO(vt_SOwaves));
        vt_soThres{cc}      = single(nm_thr);

        clear st_Cnf

        %% Detect spindle events

        fprintf('\n	** Processing fast spindles: ')
        st_Cnf              = struct;
        st_Cnf.fsampling	= nm_fsample;
        st_Cnf.minnumosc	= nm_SPINDLENUMOSC;
        st_Cnf.timebounds	= vt_timeSpindles;
        st_Cnf.hypnogram	= vt_hypnogram;
        st_Cnf.stage        = vt_stageSpindles;
        st_Cnf.rms          = vt_rmsFS;
        st_Cnf.rawEEG       = vt_signalEEG;
        st_Cnf.freqband     = vt_freqSpindleFast;
        st_Cnf.method       = 'fixed';

        mx_spindleFast	= fn_detectsleepSpindles(vt_signalFS,st_Cnf);
        clear st_Cnf

        fprintf('\n	** Processing slow spindles: ')
        st_Cnf              = struct;
        st_Cnf.fsampling	= nm_fsample;
        st_Cnf.minnumosc	= nm_SPINDLENUMOSC;
        st_Cnf.timebounds	= vt_timeSpindles;
        st_Cnf.hypnogram	= vt_hypnogram;
        st_Cnf.stage        = vt_stageSpindles;
        st_Cnf.rms          = vt_rmsSS;
        st_Cnf.rawEEG       = vt_signalEEG;
        st_Cnf.freqband     = vt_freqSpindleSlow;
        st_Cnf.window       = nm_SMOOTHWINDOW;
        st_Cnf.method       = 'fixed';

        mx_spindleSlow	= fn_detectsleepSpindles(vt_signalSS,st_Cnf);

        clear st_Cnf

        mx_spindleDet{cc,1}	= mx_spindleFast;
        mx_spindleDet{cc,2}	= mx_spindleSlow;

        %% extracting spindle feature
        if spindle ==1
            trl_mx = new_header{ff,1}.stim; trl_lim_s = nm_fsample.* [-3.04 3.04]; detected_Fspdl=[];  trl_mx2=[];
            trl_mx = trl_mx(trl_mx(:,2)==1 , 1);
            for i=1:length(trl_mx), trl_mx2(i,:) = trl_mx(i) + [trl_lim_s(1):trl_lim_s(2)]; end
            trl_mx2 = trl_mx2'; %timextrls
            trl_time = -3.04:1/nm_fsample:3.04;
            for i=1:size(mx_spindleFast,1)
                detected_Fspdl = [detected_Fspdl mx_spindleFast(i,1):mx_spindleFast(i,2)];
            end
            trl_V = trl_mx2(:); trl_V = trl_V'; spdl_pos = zeros(1,length(trl_V));
            spdl_pos(detected_Fspdl)=1;
            spdl_pos = reshape(spdl_pos, [size(trl_mx2,1),size(trl_mx2,2)]);
            spdl_pos = spdl_pos';
            % extracting spindle presence
            ft_limits = nearest(trl_time,-nm_STIMDISTANCE):nearest(trl_time,0); % extracting the spindle features in that range
            masked_spindles = spdl_pos;
            masked_spindles = masked_spindles(:,ft_limits);
            for i=1:size(masked_spindles,1)
                spindle_features{cci , ff}(i,:) = double(any(masked_spindles(i,:))); % has_spindle, yes: 1, no: 0
            end
            continue;
        end
        %% Extract time trials
        for nm_markType = 1:2
            fprintf('\n	  > Extract time trials for %i-click: ',nm_markType)
            st_Cnf              = struct;
            st_Cnf.markerTime	= mx_markTime(:,nm_markType);
            st_Cnf.fsampling    = nm_fsample;
            st_Cnf.window       = vt_timeLim;
            st_Cnf.timetrial    = vt_timeEpoch(1:2:end);

            vt_timeTrialsEEG	= fn_extracttimetrials(vt_signalEEG,st_Cnf);
            vt_timeTrialsSO     = fn_extracttimetrials(vt_signalSO,st_Cnf);
            vt_timeTrialsFS     = fn_extracttimetrials(vt_signalFS,st_Cnf);
            vt_timeTrialsSS     = fn_extracttimetrials(vt_signalSS,st_Cnf);
            vt_timeTrialsFSrms  = fn_extracttimetrials(vt_rmsFS,st_Cnf);
            vt_timeTrialsSSrms	= fn_extracttimetrials(vt_rmsSS,st_Cnf);

            switch nm_markType
                case 1
                    mx_dataClick_1{cc,1}	= mx_markTime(:,nm_markType);
                    mx_dataClick_1{cc,2}    = vt_timeTrialsEEG;
                    mx_dataClick_1{cc,3}    = vt_timeTrialsSO;
                    mx_dataClick_1{cc,4}    = vt_timeTrialsFS;
                    mx_dataClick_1{cc,5}    = vt_timeTrialsSS;
                    mx_dataClick_1{cc,6}    = vt_timeTrialsFSrms;
                    mx_dataClick_1{cc,7}    = vt_timeTrialsSSrms;
                case 2
                    mx_dataClick_2{cc,1}    = mx_markTime(:,nm_markType);
                    mx_dataClick_2{cc,2}    = vt_timeTrialsEEG;
                    mx_dataClick_2{cc,3}    = vt_timeTrialsSO;
                    mx_dataClick_2{cc,4}    = vt_timeTrialsFS;
                    mx_dataClick_2{cc,5}    = vt_timeTrialsSS;
                    mx_dataClick_2{cc,6}    = vt_timeTrialsFSrms;
                    mx_dataClick_2{cc,7}    = vt_timeTrialsSSrms;
            end
        end
        vt_timeTrial	= vt_timeEpoch(1:2:end);
        %% Compute event features
        for nm_markType = 1:2
            fprintf('\n	  > Compute event features for %i-click: ',nm_markType)
            st_Cnf              = struct;
            st_Cnf.SOlocations	= vt_SOwaves;
            st_Cnf.markerTime	= mx_markTime;
            st_Cnf.currentMark	= nm_markType;
            st_Cnf.threshold    = nm_PEAKTHRESHOLD;
            st_Cnf.stimDistance = round(nm_STIMDISTANCE * nm_fsample);
            st_Cnf.peakDistance	= round(nm_STIMDISTANCE * nm_fsample);
            st_Cnf.fSampling    = nm_fsample;

            st_Cnf.FSlocations	= mx_spindleFast;
            st_Cnf.SSlocations	= mx_spindleSlow;
            st_Cnf.FSFilter     = vt_signalFS;
            st_Cnf.SSFilter     = vt_signalSS;
            st_Cnf.FSrms        = vt_rmsFS;
            st_Cnf.SSrms        = vt_rmsSS;
            st_Cnf.hilbDelay	= nm_hilDelay;
            st_Cnf.hilbFilter	= vt_hilbFilter;
            st_Cnf.hilbPhase	= vt_phaseSO;
            st_Cnf.rawEEG       = vt_signalEEG;
            st_Cnf.spindFreq	= vt_freqSpindleFast;
            st_Cnf.timeSpind	= [0.25,1.75];
            [st_Features,st_Logs]	= fn_obtaintrialfeatures(...
                vt_signalSO,st_Cnf);
            clear st_Cnf
            switch nm_markType
                case 1
                    mx_dataClick_1{cc,10}	= st_Features;
                    mx_dataClick_1{cc,11}	= st_Logs;
                case 2
                    mx_dataClick_2{cc,10}	= st_Features;
                    mx_dataClick_2{cc,11}	= st_Logs;
            end
        end

    end 
    if spindle ==0 
        SO_features{ff,1} = mx_dataClick_1;
    end
end
if spindle ==1
    varargout{1} = spindle_features;
else
    varargout{1} = SO_features;
end

end



