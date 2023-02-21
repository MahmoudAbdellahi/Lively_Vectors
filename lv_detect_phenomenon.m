function [result] = lv_detect_phenomenon(cfg)
% takes continuous data or 3d and gets the phenomena: detects them and
% marks their timing and extracts the features of each phenomenon and this
% is determined according to the excel file so if you want to add
% phenomenon add its name and its criteria and the function should use each
% criterea and detect the event accordingly and returns a matrix filled
% with events and their features values...

% works on fieldtrip structure and the data of interest should be inside
% cfg.data.trial
%  example
% cfg.phenomenon = 'SO';
% cfg.fun_seq = [{'band'} {'duration_zcrossings'} {'spec_thresholded'}];
% cfg.data.fs = 200;

% any function could be the end on the sequence of functions and the result
% is stored in result.measures

data = cfg.data;


% 3d is now 2d because we work on one channel and squeeze before calling
% this function. So, from trl_time to vector .. and here we work on vector
if isfield(data,'trial')
    if size(data.trial{1},1)~=1 % if 1 then it's 1_time
        data.trial = reshape( permute(data.trial, [2 1]), 1,[]); %1_trltime
        data.time = cell([0 (1:length(data.trial))/data.fsample]);
    else
        data.trial = cell2mat(data.trial);
    end
end

if any(isnan(data.trial)), data.time{1,1}(isnan(data.trial))=[]; data.trial(isnan(data.trial))=[]; 
    warning('lv: data has nans and they were removed'); end

table_specs = readtable('D:\sul''s code\Matt\sleep\erps\Organised\New exp\lv_phenomenon_specs.xlsx');
titles = table_specs.Properties.VariableNames;
variables = table_specs.Variables;

% cfg.phenomenon: name of phenomenon (ex: 'SO'), specs will get the row of our phenomenon
specs = variables( ismember(variables(:,ismember(titles,'name')), cfg.phenomenon) ,:);

fprintf(['\n Detecting phenomenon: ' char(specs(1)) ' \n']);

fun_seq = cfg.fun_seq; % sequence of function execution (example: [{'band'} {'duration_zcrossings'} {'amplitude'}])
result = data; % to have the original data as initial result and  then keep updating the result squentially in the loop

if isfield(cfg,'specs')==1,specs=cfg.specs; end  % override values if we want to change from those in the sheet

specs = specs(ismember(titles,fun_seq)); % to get only the required specs


for i=1:length(fun_seq) % dynamic variable definition, starting after first one because the first is the name
    var = char(specs(i)); % the value of the variable from the excelsheet
    if ~isempty(var) % to check band is set
        eval(['result = lv_calculate_' char(fun_seq(i)) '(result,var);']);
    end
end


% TOTEST: handling the case in which we have trials so that we put the idx of the
% event in the respective trial in it's correct time (result.trial_ref_measures)
if isfield(data,'trial')
    if size(cfg.data.trial{1},1)~=1 % if the original cfg.data was 3d we will put the id of the event in the time of the trial to be able to get itt..
        if isfield(result,'measuresNames') % if we have the measure that contains flags of positions of the phenom. we will put these flags in the same format as the original data
            if strcmp(result.measuresNames(1),'index')
                id = 1:length(data.trial); % the 1d vector
                event_pos = zeros(size(data.trial));
                for i=1:size(result.measures,1) % intersect the idx of the event and the vector indices, to put the idx of event in the correct idx
                    common_pos=intersect(id, result.measures(i,1));
                    if ~isempty(common_pos), event_pos(common_pos)=i; end
                end
                result.trial_ref_measures = permute(reshape(event_pos, size(cfg.data.trial,2),[]), [2 1] );
                % result.measures{1} = permute(  reshape( result.measures{1} , size(data.trial,2) ,[]),   [2 1]);
            end
        end
    end
end
% events that are cut by the trial end won't be removed .. even if the
% trials were not continuous and there is edge artifact this artifact won't
% simulate the behavior of the phenom. for the end of one trial and the
% beg. of the other .. but what's likely is that they contained the phenom.
% when they were continuous before being cut into trials .. so we should
% leave these phenom. it's important to leave them to be able to detect a phenom happening
% at stim onset.

end

%% function for every measure ...
%%
function result = lv_calculate_band(result,var) % see the end of file for filter visualization
% vector band-pass filtering but result argument is a whole struct to be able to use ft_preprocessing
% two-pass filter by defauls (forward and backward filtering)
% fir is better in offline analysis gives better filter estimation, order in fir
% is the no. of points in time of the filter and they should be enough to
% multiply with the signal and get the similarity between them so aat least
% 3 cycles of the lowest freq.
cfg_preprocessing = [];
cfg_preprocessing.demean = 'yes';
cfg_preprocessing.bpfilter = 'yes';
cfg_preprocessing.bpfilttype = 'fir';
cfg_preprocessing.bpfreq = str2num(var);
cycles = 3; order = (cycles/min(cfg_preprocessing.bpfreq))*result.fsample;
cfg_preprocessing.bpfiltord = round(order);

result.trial = mat2cell(result.trial, 1, length(result.trial));% to cell so that ft accepts it

result= ft_preprocessing(cfg_preprocessing, result);

result.trial = cell2mat(result.trial);

result.measures = result.trial; % incase it was the last function in the seq. always look at result.measures

fprintf(['\n Done filtering in range: [' num2str(cfg_preprocessing.bpfreq) '] \n']);
end

%%
function result = lv_calculate_duration_zcrossings(result,var)
% gets the positive to negative zero-crossings and it works on result.trial


duration_samples = round(result.fsample* str2num(var) );

pos_to_neg_zcrossings = (diff(sign(result.trial)) <= -1);

event_sample_marker = find(pos_to_neg_zcrossings==1);

consecutive_crossings_samples = diff( event_sample_marker ); % no. samples between consecutive crossings

if any( consecutive_crossings_samples ==1 ) % 1 0 0 -1 meaning to zeros were found that's a problem because
    % the signal should not shift like that into negative we would normally expect 1 0 -1 or 1 -1
    error('lv: found a signal with two consecutive points at zero when it was transitioning from positive to negative!');
end


good_event = find( consecutive_crossings_samples>=duration_samples(1) &  consecutive_crossings_samples<=duration_samples(2));

for i=1:length(good_event)
    % first crossing sample and the consecutive is the second crossing
    idx{i,1} = [event_sample_marker(good_event(i)): event_sample_marker(good_event(i)+1)]; % good events
end

temp = zeros(size(result.trial)); temp( cat(2,idx{:}) )=1;

result.trial = result.trial.*temp; % the output of this one is the masked signal according to the phenom.
result.measures = result.trial;
result.idx = idx; % very important, this contains the indices of all phenomena

fprintf(['\n Done thresholding on zero_crossings. \n']);
% visualize phenomenon here
% result.trial = mat2cell(result.trial, 1, length(result.trial));
% artf=ft_databrowser([],result);
end


%%
function result = lv_calculate_spec_thresholded(result,var)
% calculates the measure like max to min amplitude or moving power and does thresholding and it works on result.trial

var = split(var);

% the event should be determined before (if like SO) going here because how will you
% calculate amplitude to unknown phenom.

switch var{1} % always the name of the method
    case 'max_min'
        % getting measures from every event
        for i=1:length(result.idx)
            event_signal{i,1} = result.trial( result.idx{i} );
        end
        amplitudes = ( cellfun(@(x) (max(x)-min(x)),event_signal,'Un',0) );
        max_vals = ( cellfun(@(x) (max(x)),event_signal,'Un',0) );
        min_vals = ( cellfun(@(x) (min(x)),event_signal,'Un',0) );
        
        if length(var)==3 % percentile
            if strcmp(var{3},'p')==1
                threshold=prctile(cell2mat(amplitudes) ,str2num(var{2})); end
        end
        % to accept all SOs put threshold=0
        threshold=0;
        good_event = cell2mat(amplitudes)>=threshold;
        result.measures = [result.idx(good_event) event_signal(good_event) ... % index and the signal(filtered in phenom band) itself stored, so you can timelock later on let's say trough you know the idx, make it event and define trls
            amplitudes(good_event) max_vals(good_event) min_vals(good_event)]; % the amplitude, max, min .. if you want you can another stat here
        % just execute it in the cellfun and add it here, or even later because we have the signal
        result.measuresNames = [{'index'} {'event_signal'} {'max_min amplitudes'} {'max'} {'min'}];
        
        % visualize phenomenon here, beaware that to visualize we change
        % result.trial so don't look at result.trial in the return value
        %         temp = zeros(size(result.trial));
        %         temp(cell2mat(result.idx)) = result.trial(cell2mat(result.idx));
        %         result.trial = mat2cell(temp, 1, length(result.trial));
        %         artf=ft_databrowser([],result);
        
        % to make sure to see the time_locked SO to be sure that there is
        % no filtering induced time shifts between the raw and filtered data
        warning('lv: please see the trough locked signals before proceeding with SO detection if there is time shift, consider changing filter band/order');
        
    case {'rms','hilbert'} % rms is different because we will need the continuous signal to estimate power using a sliding window
        % example in sheet: rms 0.2 75 p 0.5 3, method, moving wind,
        % no. percentile, percentile, min duration, max duration
        if strcmp(var{1},'rms')==1
            win_time = str2num(var{2}); win_s = result.fsample*win_time;
            ws = floor(win_s/2);
            rms_pw = zeros(size(result.trial));
            % getting windows 1 sample shifts without looping
            id = (1:length(result.trial)-ws*2)';
            id = [id id+(1:(ws*2))];
            temp = result.trial(id);
            
            % rms calculation on every window
            pw = (sqrt(mean((temp.^2), 2)))'; % this variable to use later on when getting the percentile because percentile stats. shouldn't include zeros
            rms_pw(ws+1 : length(result.trial)-ws) = pw;
            clear temp id;
            extracted_pw = rms_pw;
        end
        if strcmp(var{1},'hilbert')==1
            hilbert_dat = hilbert( result.trial(:) ); % time in first dimension
            [pw,extracted_pw] = deal( (abs(hilbert_dat).^2)' );
        end
        
        if length(var)==6 % percentile
            if strcmp(var{4},'p')==1
                threshold=prctile(pw ,str2num(var{3})); end
        end
        
        pos_good = extracted_pw > threshold;
        
        CC = bwconncomp(pos_good);
        
        win_time = [str2num(var{5}) str2num(var{6})]; win_s = result.fsample*win_time;
        event = cell2mat(cellfun(@(x) (length(x)>win_s(1) & length(x)<win_s(2)), CC.PixelIdxList, 'Un',0));
        result.idx = CC.PixelIdxList(event)';
        
        % getting measures from every event
        for i=1:length(result.idx)
            event_signal{i,1} = result.trial( result.idx{i} );
            event_rms_pw{i,1} = extracted_pw( result.idx{i} );
        end
        max_vals = ( cellfun(@(x) (max(x)),event_signal,'Un',0) );
        min_vals = ( cellfun(@(x) (min(x)),event_signal,'Un',0) );
        
        
        result.measures = [result.idx event_signal ... % index and the signal(filtered in phenom band) itself stored, so you can timelock later on let's say trough you know the idx, make it event and define trls
            event_rms_pw max_vals min_vals]; % the rms_power, max, min .. if you want you can another stat here
        % just execute it in the cellfun and add it here, or even later because we have the signal
        result.measuresNames = [{'index'} {'event_signal'} {'event_rms_pw'} {'max'} {'min'}];
        
        % visualize phenomenon here
        %         temp = zeros(size(result.trial));
        %         temp(cell2mat(result.idx)) = result.trial(cell2mat(result.idx));
        %         result.trial = mat2cell(temp, 1, length(result.trial));
        %         artf=ft_databrowser([],result);
        
        
    case 'instant_phase_magnitude' % instantenuous magnitude and phase using hilbert transform
        hilb_analytical = hilbert(result.trial(:)); hilb_analytical=hilb_analytical';
        result.measures = [mat2cell(1:length(result.trial),1,length(result.trial)) mat2cell(result.trial,1,length(result.trial)) ...
            mat2cell(abs(hilb_analytical),1,length(result.trial))  mat2cell(angle(hilb_analytical),1,length(result.trial))];
        result.measuresNames = [{'index'} {'event_signal'} {'magnitude'} {'phase'}]; % here we have 'event' at every time point because we just assume that this is an event
        % to get the phase locked to stim here we have the phase of the
        % original signal with the high sampling rate you can still see the
        % samples and oversample the time vector of the trial and get the
        % circmean of the phase values in the time durations you want
end

fprintf(['\n Done applying specification thresholding and extracting phenomenon measures using: ' var{1} ' \n']);




end


%%
function result = lv_calculate_coupled_SOspindle(result,var) 

if strcmp(var,'modulation_idx')==1 % this is actual phase amplitude coupling
    % it can be done between any two phenom .. read lv_examples for details
    % and example
    % the point is that coupling strength isn't real until we subtract the
    % surrogate data from it because the consistency is always existing and
    % always positive so the real consistency is the consistency vs
    % surrogate so it's like subtraccting chance level .. and because we are using the magnitude of the faster
    % phenom. so the higher the magnitude the higher the contribution to
    % the mean so the coupling strength that's why we say that the phase is
    % the phase of the slower phenom with highest amplitude of the faster
    % phenom because the highest amplitudes get to contribute more to the
    % mean 
    
    % 5od hena el so measures w spindle measures w el hilbert bta3hom elle
    % tele3 w hena (kol el calls fe examples) ttba2 kalam el coupling ....
    % netala3 el indices bta3t el phenom w nshof el hilbert info el
    % corresponding w b3d kda ne measure el coupling 3la 7asab el below
    % code
    so_measures = result.so_measures;
    spindle_measures = result.spindle_measures;
    so_hilbert = result.so_hilbert;
    spindle_hilbert = result.spindle_hilbert;
    
%     we need another intersection here with the correct trials continuous
%     time to be able to get the exact phase of SO_spindle_correct and
%     SO_spindle_incorrect .. and compare the number of these individual
%     events for correct vs. incorrect.. if you want the numebr not
%     modulation take the intersection in the main code not in this
%     function.. you can normalise the events by dividing by the total
%     SO_spindle_correct events despite the phase transition to make it
%     probability and also may be able to compare it to 0.5 chance level..
%     the event may be good because the reactivation happens on the trough
%     of spindles so the events (timepts wise) could be a good indicator vs
%     taking the whole trial because within one trial it might be correct
%     (trough of spindles) then random (peak of spindles) but taking the
%     correctness instersected with SO_spindles at every time pt should be
%     good

    events = cell2mat(spindle_measures.measures(:,1)); % events/pts to look at 
%     for i=1:size(spindle_measures.measures,1), [~,idx] = max(spindle_measures.measures{i,2}); events(i)=spindle_measures.measures{i,1}(idx); end
  
    so_spindle_idx = intersect(cell2mat(so_measures.measures(:,1)'), events); 
    
    if isfield(result,'idx_limiter') % this idx_limiter is an additional condition on the extracted events 
        % so it can be the correct trials idx in time so that we test the
        % PAC for these particular regions
        so_spindle_idx = intersect(so_spindle_idx, result.idx_limiter);
    end
    
    amplitude = spindle_hilbert.measures{1, 3}(so_spindle_idx); % magnitude of the faster phenom.
    phase = so_hilbert.measures{1, 4}(so_spindle_idx); % phase of the slower phenom.
    
    amplitude = reshape(amplitude,1,[]); phase = reshape(phase,1,[]);
    % code from Canolty 2006
    srate=0; %% sampling rate used in this study, in Hz .. you may put it 0 if the number of events is really low because you won't like to cut any data !
    numpoints=length(amplitude); %% number of sample points in raw signal
    numsurrogate=10000; %% number of surrogate values to compare to actual value
    minskip=srate; %% time lag must be at least this big
    maxskip=numpoints-srate; %% time lag must be smaller than this
    skip=ceil(numpoints.*rand(numsurrogate*2,1));
    skip(find(skip>maxskip))=[];
    skip(find(skip<minskip))=[];
    skip=skip(1:numsurrogate,1);
    surrogate_m=zeros(numsurrogate,1);
    % HG analytic amplitude time series, uses eegfilt.m from EEGLAB toolbox
    %amplitude=abs(hilbert(eegfilt(x,srate,80,150))); % mo because we already have our magnitude and phase 
    % theta analytic phase time series, uses EEGLAB toolbox
    %phase=angle(hilbert(eegfilt(x,srate,4,8))); 
    % complex-valued composite signal
    z=amplitude.*exp(1i*phase);
    % mean of z over time, prenormalized value
    m_raw=mean(z);
    % compute surrogate values
    for s=1:numsurrogate
        surrogate_amplitude=[amplitude(skip(s):end) amplitude(1:skip(s)-1)];
        surrogate_m(s)=abs(mean(surrogate_amplitude.*exp(1i*phase)));
    end
    % fit gaussian to surrogate data, uses normfit.m from MATLAB Statistics toolbox
    [surrogate_mean,surrogate_std]=normfit(surrogate_m);
    % normalize length using surrogate data (z-score)
    m_norm_length=(abs(m_raw)-surrogate_mean)/surrogate_std;
    m_norm_phase=angle(m_raw);
    m_norm=m_norm_length*exp(1i*m_norm_phase); % modulation index which is complex the magnitude of it is the 
    % coupling strength and the phase is the coupling phase
    
    result.coupling_strength = m_norm_length; % zvalue that can be converted to pval
    result.phase_of_coupling = m_norm_phase; % the (slower) phase with the highest (faster) magnitude values  
    % not sure that this is the actual coupling phase because the coupling
    % phase analysis is different and should consider the magnitude like
    % taking the circular mean of amplitudes or binning the phase angles or transitions

    % will keep the coupling uncorrected because that's the actual coupling value 
    % and report the pvalue alone not as zvalue because zvalue is not
    % coupling but more of a measure of significance
    result.coupling_strength = abs(m_raw);
    result.pval = mean(surrogate_m>abs(m_raw));

elseif strcmp(var,'max_spindle')==1 
    % this method is the one in Endogenous memory reactivation during
    % sleep.. and it is based on getting the max of every spindle and then
    % looking at the phase angle of that and then checking if this will
    % make a non uniform distribution because if so then it will mean that
    % the angles are concentrated in specific phase..
    so_measures = result.so_measures;
    spindle_measures = result.spindle_measures;
    so_hilbert = result.so_hilbert;
    spindle_hilbert = result.spindle_hilbert;

    [amplitude,id] = cellfun(@(x) max(spindle_hilbert.measures{1, 3}(x)), spindle_measures.measures(:,1));  % max of spindle
    for i=1:length(id), ids(i) = spindle_measures.measures{i,1}(id(i)); end % the actual id on the continuous vector
    % SO_spindle event
    SO_spindle = intersect(ids, cell2mat(so_measures.measures(:,1)')); 
    SO_spindle = cell2mat(so_measures.measures(:,1)'); % for the phase analysis of correct pts 

    phase = so_hilbert.measures{1, 4}(SO_spindle); % phase of the slower phenom.
    %circ_plot(phase','pretty','bo',true,'linewidth',2,'color','b')
    [pval z] = circ_rtest(phase);
    result.uniformity = pval; 
    result.phase = circ_mean(phase');

    result.phase = rad2deg(phase); result.phase(result.phase<0)= result.phase(result.phase<0)+360; % for phase count
    temp = result.phase;
    result.phase = [mean(temp>180 & temp<360)  mean(temp>0 & temp<180)]; % up down
%     result.phase = [mean(temp>0&temp<90 | temp>270&temp<360)  mean(temp>90&temp<270)]; % positive negative -halves

elseif strcmp(var,'fieldtrip')==1 % that's fieldtrip's implementation for canolty mean vector
    % length (mvl) and Tort's modulation index based on histogram bins 

    error('check page 205 of Notes 9'); % because this detects the phenom as well so we shouldn't use it after pehnom detection
    error('there is another one in feature_extractor')
    so_measures = result.so_measures;
    spindle_measures = result.spindle_measures;
    so_hilbert = result.so_hilbert;
    spindle_hilbert = result.spindle_hilbert; 

    events = cell2mat(spindle_measures.measures(:,1)); % events/pts to look at 

    so_spindle_idx = intersect(cell2mat(so_measures.measures(:,1)'), events); 
    
    amplitude = (spindle_hilbert.measures{1, 3}(so_spindle_idx))'; % magnitude of the faster phenom.
    phase = so_hilbert.measures{1, 4}(so_spindle_idx); % phase of the slower phenom.
     
    % crossfreq = ft_crossfrequencyanalysis(cfg, freqlo, freqhi); % freqlow is the slower phenom which is SO
    % calling the function will require data formatting so instead we will
    % get the low level function from the cossfrequency function and use
    % them here .. and get other low lvl function if other methods are needed
    
    % phase/amplitude is a matrix with each row representing a freq. phase/amplitude
    option = 2; % 1 for Canolty's method and 2 for Tort's
    
    if option==1
        result.coupling_strength = data2mvl(phase,amplitude); % this reshape is in case of one row
    elseif option==2
        nbin= 20;
        [result.coupling_strength result.phase_of_coupling] = data2pac(phase,amplitude, nbin); % 20 is the number of bins (in tort's method they set it to 18 and ft to 20) 
    end
    % Permutation to get the pvalue of coupling strength and know if it's significant
    % the permuted amplitude-phase assignment is done via cutting the amplitude in half at random pts. and re-run
    permutations = 1000; rng(1);
    sp = randi([2,size(phase,2)-1], permutations,1); % random splits for amplitude to create new amplitude-phase data
    for i=1:size(sp,1)
        amplitude_temp = [amplitude(:,sp(i,1):end) amplitude(:,1:sp(i,1)-1)];
        if option==1
            coupling_strength_temp(:,:,i) = data2mvl(phase,amplitude_temp); % this reshape is in case of one row
        elseif option==2 
            coupling_strength_temp(:,:,i) = data2pac(phase,amplitude_temp, nbin); % 20 is the number of bins in tort's method
        end
    end
    % Pval
    for i=1:size(result.coupling_strength,1)
        for j=1:size(result.coupling_strength,2)
            result.pval(i,j) = mean(squeeze(coupling_strength_temp(i,j,:)) > result.coupling_strength(i,j));
        end 
    end
else
    % not actual PAC coupling it's cooccuring SO and spindles (so spindle complexes) takes the timing of coupling and w
    % here we need the id of SO and spindles to calculate when they are coupled
    % and return that as a new event..
    duration = str2num(var);
    duration_samples = duration*result.so_measures.fsample;
    
    [~,mini] = cellfun(@(x) min(x), result.so_measures.measures(:,2),'Un',0);
    min_id = cell2mat(cellfun(@(x,y) y(x), mini,result.so_measures.measures(:,1),'Un',0)); %  to see the real id in samples
    
    for i=1:length(min_id), search_proximity(i,:) = [min_id(i) : min_id(i)+duration_samples]; end
    
    spindles_ids = result.spindle_measures.measures(:,1); spindles_ids=cell2mat(spindles_ids);
    event=[];
    for i=1:size(search_proximity,1) % searching for co-occurring spindles within SO trough proximity
        if ~isempty( intersect(spindles_ids, search_proximity(i,:)) )
            event=[event i];
        end
    end
    
    result.so_spindle_complex.measures=result.so_measures.measures(event,:); % SOs that have spindles
    
    fprintf(['\n Done extracting SO_spindles complexes by searching spindles after ' var ' from the SO trough \n']);
    
    
end

end


% %% visualize filter
% % get filter's b from the fil1 function or other functions that create the
% % filter b is the filter in time so we want to transform it to frequency
% % and plot it to see how good it fits our desired frequencies.. code from
% % mike lec:125
% plot(B); % just the kernel 'b' in time
% filtpow = abs(fft(B)).^2;
% nyq = 200/2;
% hz = linspace(0,nyq,floor(length(filtpow)/2)+1); % to half sampling and the +1 for 0hz
% filtpow = filtpow(1:length(hz)); % cutting to freq of interest
% plot(hz(1:22), filtpow(1:22)) % to see the band of interest more clearly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION from fieldtrip's crossfrequency for PAC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mvldata] = data2mvl(LFsigtemp,HFsigtemp) % Canolty 2006
% calculate  mean vector length (complex value) per trial
% mvldata dim: LF*HF
% each row is a different frequency for phase inside LF for which you want
% PAC and also each row is a different frequency for amplitude inside HF
LFphas   = LFsigtemp;
HFamp    = HFsigtemp;
mvldata  = zeros(size(LFsigtemp,1),size(HFsigtemp,1));    % mean vector length

for i = 1:size(LFsigtemp,1)
  for j = 1:size(HFsigtemp,1)
    mvldata(i,j) = nanmean(HFamp(j,:).*exp(1i*LFphas(i,:)));
  end
end

mvldata = abs(mvldata); % the magnitude of the modulation index is the coupling strength

end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pacdata,lv_phase_angle] = data2pac(LFsigtemp,HFsigtemp,nbin)
% calculate phase amplitude distribution per trial
% pacdata dim: LF*HF*Phasebin

pacdata = zeros(size(LFsigtemp,1),size(HFsigtemp,1),nbin);

Ang  = LFsigtemp;
Amp  = HFsigtemp;
[dum,bin] = histc(Ang, linspace(-pi,pi,nbin));  % binned low frequency phase
binamp = zeros (size(HFsigtemp,1),nbin);      % binned amplitude

for i = 1:size(Ang,1)
  for k = 1:nbin
    idx = (bin(i,:)==k);
    binamp(:,k) = mean(Amp(:,idx),2);
  end
  pacdata(i,:,:) = binamp;
end

% lv: getting the phase of the amplitudes
bins = linspace(-pi,pi,nbin)';
% we convert the amplitudes into weights for the binned angles so that the circular mean of their amplitude weighted values is the phase of coupling
w = squeeze(pacdata)./nansum(squeeze(pacdata));
bins(isnan(w))=[]; w(isnan(w))=[];
lv_phase_angle = circ_mean(bins, w);


% now the KL distances from the uniform distribution for every point
nlf = size(Ang, 1); nhf = size(Amp, 1);
% from fieldtrip, the part of getting the KL distance
Q =ones(nbin,1)/nbin; % uniform distribution
mi = zeros(nlf,nhf);
for i=1:nlf
    for j=1:nhf
        P = squeeze(pacdata(i,j,:))/ nansum(pacdata(i,j,:));  % normalized distribution
        mi(i,j) = nansum(P.* log2(P./Q))./log2(nbin); % KL distance
    end
end
pacdata = mi;

end % function
