function [ good_trials ] = lv_reduce_trials(cfg)
% good_trials are always sorted from good to bad according to method..
% reduces the trials based on measures .. stats etc.,
% cfg includes the data and method and could include many datasets
% according to the method used ..

% As a matter of fact this function can be used to reduce any dimension not
% just for trials ...

% returing idx of good trials to be able to use it for sorting trials and
% color coded erps and visualization as well ..

options = {'var','rms power','mean_centroid','low var',...
    'low rms power'};

if ~isfield(cfg,'method')
    idx = listdlg('ListString',options, 'ListSize',[450,450]);
    cfg.method = char( options(idx) );
end

% fprintf(['\n Reducing features based on: ' cfg.method ' \n']);


data = cfg.data;
switch cfg.method
    case {'var','rms power'} % these measures compress time and compare to noise .. so data vs noise and get the data far from noise
        noise = cfg.noise; % we assume a second data set that contains other trials that can be
        % the same data set or wake or control night etc., notably it is
        % used to get the cenroid of its trials to be noise.
        if strcmp(cfg.method,'var'), noise_measure = squeeze(var(noise.trial,0, 3)); data_measure = squeeze(var(data.trial,0, 3)); end % trl_ch
        if strcmp(cfg.method,'rms power'), noise_measure = squeeze(mean(noise.trial.^2, 3)); data_measure = squeeze(mean(data.trial.^2, 3)); end % trl_ch

        noise_centroid = squeeze(mean(noise_measure,1)); % 1_ch
        distance = pdist2(noise_centroid, data_measure,'euclidean'); % acts on rows and every pair.. can use spearman,, performs on all trials at once
        [~, good_trials] = sort( distance,'descend' ); % good are the ones away from noise

    case 'mean_centroid' % this doesn't compress time
        noise = cfg.noise; % we assume a second data set that contains other trials that can be
        % the same data set or wake or control night etc., notably it is
        % used to get the centroid of its trials to be noise. .. if noise is
        % one channel repeat it before calling the function so channels
        % will match ... manipulate .trial if you want away from flat etc.,
        noise_centroid = squeeze(mean(noise.trial,1)); % ch_time
        %noise_centroid  = median(noise_centroid ,1); %if you want to compress on channels as well
        for ch=1:size(noise_centroid,1)
            distance_temp(ch,:) = pdist2(noise_centroid(ch,:), squeeze(data.trial(:,ch,:)),'euclidean'); % acts on rows and every pair.. can use spearman,, performs on all trials at once
        end
        [~, good_trials] = sort( median(distance_temp,1),'descend' ); % good are the ones away from noise

        % choosing the no. trials using iqr
        %iqr_val = iqr( median(distance_temp,1)); q1=prctile( median(distance_temp,1),25);
        %good_trials = find( median(distance_temp,1) > q1-(iqr_val));

    case {'low var','low rms power'} % these measures compress time and consider low var and low power as bad trials
        if strcmp(cfg.method,'low var'), trls_measure = median( var(data.trial, 0,3), 2); end % var across time and median of channels
        if strcmp(cfg.method,'low rms power'), trls_measure = median( mean(data.trial.^2, 3), 2); end % var across time and median of channels

        [~,good_trials] = sort(trls_measure,'descend'); % good are with high variance/power


    case 'pre-classification' % iqr based outlier rejection
        % this one takes 2d and loops over 2nd dim and makes a adistribution with 1st dim and
        % returns the goodFlags_matrix in good_trials
        all_measure = data.trial; % ex: trl_ch
        goodFlags_matrix = nan(size(all_measure));
        for i=1:size(all_measure,2)
            % interquartile range (IQR) based outliers rejection
            iqr_val = iqr(all_measure(:,i)); q1=prctile(all_measure(:,i),25); q3=prctile(all_measure(:,i),75);
            var_threshold = [q1-(1.5*iqr_val)   q3+(1.5*iqr_val)]; % two sided threshold

            goodFlags_matrix(:,i) = all_measure(:,i)>var_threshold(1) & all_measure(:,i)<var_threshold(2);
        end
        good_trials = goodFlags_matrix;
        % trl_ch of flag values is better than other return values because
        % we will be able to see exactly where it happened and can do '&' with other measure if we want
    case 'fixing_outliers'
        % gets the outliers in a similar way as pre-classification so it
        % uses median and quartile ranges the difference is that it uses a
        % method called 'winsorization' this method sets the outliers>1.5iqr+q3 to
        % be at max value closest(not fixed threshold actually) and the ones lower q1-1.5iqr to be min value closest .. will reshape the
        % data to act on the 1st dim so we work on mdim data ...
        data=data.trial; sz = size(data);
        all_measure = reshape(data,sz(1),[]);
        for i=1:size(all_measure,2) % interquartile range (IQR) based outliers detection
            iqr_val = iqr(all_measure(:,i)); q1=prctile(all_measure(:,i),25); q3=prctile(all_measure(:,i),75);
            threshold = [q1-(1.5*iqr_val)   q3+(1.5*iqr_val)]; % two sided threshold
            lowest_val = min(all_measure( all_measure(:,i)>threshold(1) ,i )); % lowest non outlier value
            largest_val = max(all_measure( all_measure(:,i)<threshold(2) ,i )); % largest non outlier value
            all_measure(all_measure(:,i)<=threshold(1),i) = lowest_val;
            all_measure(all_measure(:,i)>=threshold(2),i) = largest_val;
            lv_progress(i,size(all_measure,2),'fixing_outliers: ');
        end
        good_trials = reshape(all_measure,sz); % this should contain the corrected trials
    case 'naive bayes based'
        % rejecting outliers that have high posterior probability of being
        % in the class that contains all trials ... we use the data from
        % noise and make the gaussian distribution and reject the trials behaving like them
        % just be careful that it's not real probability and you need a
        % threshold.. you can use it with 2d data by making the 3rd dim=1 so just ignore the 3rd dim
        noise = cfg.noise;

        feats = reshape(noise.trial, size(noise.trial,1),[]);
        data_feats = reshape(data.trial, size(data.trial,1),[]);
        mu = mean(feats,1);
        sigma = std(feats,[],1);

        for i=1:length(mu) % no prior because it's 1class we are estimating the likelihood here
            likelihood(:,i) = pdf('Normal',data_feats(:,i) ,mu(i),sigma(i));
        end

        [~,good_trials] = sort(prod(likelihood,2),'ascend'); % they are muliplied because they are happening together meaning that we get value1 for feature1 and value 2
        % for feature2 at the same time .. imagine tow events the
        % probability of each is 0.5 then the probability of them hapening
        % together is 0.25

        % because the prev. one vanishes we are assuming and trying this one
        temp = mean(squeeze(prod( reshape(likelihood, size(data.trial)) ,2)),2); % taking the mean in time not sure it will be fine i am just assuming that time pts are dependent and
        % channels are the features so they control pdf but the time pts
        % don't so we multiply probabilities of features but get the mean
        % of time
        [~,good_trials] = sort(temp,'ascend');

    case 'decorrelated'
        % loops over trials and marks a trial as good if it is decorrelated
        % between the right and left motor area.. will assume that dim2
        % contains channels and they are put in one area then the other
        % area is put in cfg.noise ... if you have a channel that can be the
        % noise then repeat it in cfg.noise to match the channels of data
        % and you will have every channel vs noise
        for i=1:size(data.trial,1)
            rho = corr(data.trial(i,:,:)',  cfg.noise.trial(i,:,:)', 'type','spearman'); % ex: 1st row the frst col from sleep corr with all trls from wake
            rhos(i,1)=mean( diag(rho) );
        end
        [~,good_trials] = sort(rhos,'ascend');


    case 'n100'
        % n100 or phenomenon
        % keep the trials with the sound erp component .. because if the
        % brain heard the sound then it should process the information of
        % the cue .. otherwise it won't,, the same thing can be applied to
        % other effects that are bounded in time not just n100 so something
        % like expected spindle after 1 sec. etc.
        % one-class classification with leave one out and verifies the
        % point in sliding window so vs all windows ..
        n100_time_bounds = [0.070 0.120]; % n100 [0.070 0.120] or any other phenomenon
        if strcmp(cfg.routine,'training')
            % training routine to build the model
            data = cfg.data; % struct array of sbjs .. each contains 3d trl_ch_time
            timeax = cfg.time; % from one sbj just the time vector
            % channel(s) chosen for phenomenon are the ones kept (the choice is given so that we don't have to deal with all ch and save memory)
            for i=1:length(data)
                idx=ones(length(data),1); idx(i)=0; idx=logical(idx);
                TRAIN = cat(1,data{idx});
                TEST = cat(1,data{~idx});

                n100_idx = nearest(timeax,n100_time_bounds(1)):nearest(timeax,n100_time_bounds(2));
                TRAIN = TRAIN(:,:, n100_idx); % cutting to n100 time and channel
                TRAIN = reshape(TRAIN, size(TRAIN,1),[]);% trl_n100time because we have one channel
                TEST = squeeze(TEST); % trl_Longtime

                fprintf(['\n lv: Building phenomenon model... sbj:' num2str(i) '\n']);
                corrModel = TRAIN; % assuming that it's model but actually they are just training trls because we are doing correlations

                % sliding on test time and apply model
                for t=1:size(TEST,2)-(length(n100_idx)-1)
                    mu_rho(i,t) = predict_corr(corrModel, TEST(:,t:t+(length(n100_idx)-1) ) ); % apply model on test time window

                    ctime_ax(1,t) = mean(timeax(t:t+(length(n100_idx)-1) ));
                    if t==n100_idx(1), n100_timept(i,1)=mu_rho(i,t);  end % to know the n100 point because it should be the highest if the classifier is working
                    lv_progress(t, size(TEST,2)-(length(n100_idx)-1), 'sliding in time> ')
                end
            end
            mo_pretty_errorbar(ctime_ax,mu_rho,(mu_rho.*0)+(median(mu_rho(:))), 0); % pt wise stats

            % if performance is good train with all and save model
            TRAIN = cat(1,data{:});
            n100_idx = nearest(timeax,n100_time_bounds(1)):nearest(timeax,n100_time_bounds(2));
            n100_corrModel = squeeze(TRAIN(:,:, n100_idx));
            save n100_corrModel n100_corrModel;

        else
            % testing routine
            load n100_corrModel n100_corrModel;
            timeax = cfg.time;
            n100_idx = nearest(timeax,n100_time_bounds(1)):nearest(timeax,n100_time_bounds(2));
            TEST = cfg.data(:,:, n100_idx); % cfg.data is now test, cutting to n100 time and channel
            TEST = reshape(TEST, size(TEST,1),[]);

            rho = corr(TEST', n100_corrModel', 'type','spearman'); % ex: 1st row the frst col from sleep corr with all trls from wake
            mu_rho = mean(rho, 2);

            [~,good_trials] = sort(mu_rho,'descend');% likelihood of having n100, higher better
        end
    case 'up_going'
        % chooses the trials with the beginning of trials falling on the
        % upgoing phase of the SO
        phenom_ch = cfg.phenom_ch;
        dat = cfg.data; sz = size(dat.trial); cfg_hold = cfg;
        dat.sampleinfo = (1:length(dat.time):length(dat.time)*sz(1))';
        dat.sampleinfo = [dat.sampleinfo dat.sampleinfo+length(dat.time)-1];
        cfg = []; cfg.channel = phenom_ch; % to channel of phenomenon
        dat = ft_selectdata(cfg,dat);
        temp = squeeze(dat.trial); dat.trial=[];
        for i=1:size(temp,1), dat.trial{1,i}=temp(i,:); end % ch_time

        % SO
        cfg=[]; data=[];
        cfg.phenomenon = 'SO';
        cfg.fun_seq = [{'band'} {'duration_zcrossings'} {'spec_thresholded'}];
        cfg.data = dat; cfg.data.fsample = 200;
        cfg.data.time = mat2cell(0:1/cfg.data.fsample:(length(cfg.data.trial)*...
            length(cfg.data.time)/cfg.data.fsample)-(1/cfg.data.fsample), 1,length(cfg.data.trial)*length(cfg.data.time));
        cfg.data=rmfield(cfg.data,[{'sampleinfo'} {'trialinfo'}]);  % we don't need them inside the function to be able to filter continuous with no issues
        [data.so_measures] = lv_detect_phenomenon(cfg);
        % marking phenom events on trials
        ids = (cell2mat(data.so_measures.measures(:,1)'))';
        SO_vec_pos = zeros(1,length(cfg.data.time{1, 1})); SO_vec_pos(ids) = 1;

        % in case we want to lock to trough
        if isfield(cfg_hold,'lock_to_trough')
            SO_trough_pos = zeros(1,length(cfg.data.time{1, 1}));
            for i=1:size(data.so_measures.measures,1), [~,id]=min(data.so_measures.measures{i,2});
                pos = data.so_measures.measures{i,1}(id);
                SO_trough_pos(pos) = 1;
            end
            good_trials = SO_trough_pos;
            good_trials = (reshape(good_trials,[sz(3) sz(1)]))'; return;
        end
        %         % visualising the trough locked signals
        %         tempo = cell2mat(dat.trial); dd = find(SO_trough_pos==1);
        %         for j=1:sum(SO_trough_pos)-1
        %             signal(j,:) = tempo( dd(j)-400 : dd(j)+400 );
        %         end
        %         plot(-400:400,mean(signal,1));


        % getting the SO phase of the found event ids
        % SO_hilbert
        cfg=[];
        cfg.phenomenon = 'Instant analytical';
        cfg.fun_seq = [{'band'}  {'spec_thresholded'}];
        cfg.specs = [{'Instant analytical'} {'0.5 2'} {''} {'instant_phase_magnitude'} {'no'}]; % override the record with new values 0.16 1.25
        cfg.data=dat; cfg.data.fsample = 200; cfg.data.time = mat2cell(0:1/cfg.data.fsample:(length(cfg.data.trial)*...
            length(cfg.data.time)/cfg.data.fsample)-(1/cfg.data.fsample), 1,length(cfg.data.trial)*length(cfg.data.time));
        cfg.data=rmfield(cfg.data,[{'sampleinfo'} {'trialinfo'}]);  % we don't need them inside the function to be able to filter continuous with no issues
        so_hilbert= lv_detect_phenomenon(cfg); % to get the phase of SO

        events_phase = rad2deg(so_hilbert.measures{1,4});
        events_phase(events_phase<0) = events_phase(events_phase<0)+360;
        if sum(events_phase>360)>0, error('strange phase angle..'); end
        up_going = [180 360];% [0 180];
        good_trials = events_phase>up_going(1) & events_phase<up_going(2);
        good_trials = good_trials & SO_vec_pos;
        good_trials=double(good_trials);
%         % downgoing code:nan
%         up_going = [0 180];
%         good_trials2 = events_phase>up_going(1) & events_phase<up_going(2);
%         good_trials2 = good_trials2 & SO_vec_pos; good_trials(good_trials2==1)=nan;


        good_trials = (reshape(good_trials,[sz(3) sz(1)]))'; % no channels this is trl_time
        % then take the first value from every trial as the mark for
        % whether TMR was on up_going phase or not
%                 time0 = nearest(dat.time  , 0);
%                 tmr_time = 0.2*200;
%                 good_trials = sum(good_trials(:,time0:time0+tmr_time),2) > round(tmr_time/2); % majority of TMR on this phase
% %                 good_trials = good_trials(:,time0); % at time 0 only
        % phase of TMR with circular mean
        temp_ph = so_hilbert.measures{1,4};
        temp_ph = (reshape(temp_ph,[sz(3) sz(1)]))';
        good_trials=circ_mean(temp_ph(:, nearest(dat.time,0):nearest(dat.time,0.2)),[], 2);
        events_phase = rad2deg(good_trials); events_phase(events_phase<0) = events_phase(events_phase<0)+360;
        if sum(events_phase>360)>0, error('strange phase angle..'); end
        good_trials=events_phase>180 & events_phase<360;

    case 'ASC'
        % Active system consolidation testing by looking at SO phase and
        % searching for spindles that comes after this SO because those
        % thalamo-cortical spindles should carry ripples and reactivation
        % on their troughs .. so we return trl_time of flags showing the
        % beginning of that spindle as 1 .. so you can then limit the
        % trials to those with that has this spindle or time lock to the this time
        % pt and be spindle locked to the beginning of spindles ..
        phenom_ch = cfg.phenom_ch;
        dat = cfg.data; sz = size(dat.trial); cfg_hold = cfg;
        dat.sampleinfo = (1:length(dat.time):length(dat.time)*sz(1))';
        dat.sampleinfo = [dat.sampleinfo dat.sampleinfo+length(dat.time)-1];
        cfg = []; cfg.channel = phenom_ch; % to channel of phenomenon
        dat = ft_selectdata(cfg,dat);
        temp = squeeze(dat.trial); dat.trial=[];
        for i=1:size(temp,1), dat.trial{1,i}=temp(i,:); end % ch_time

        % SO
        cfg=[]; data=[];
        cfg.phenomenon = 'SO';
        cfg.fun_seq = [{'band'} {'duration_zcrossings'} {'spec_thresholded'}];
        cfg.data = dat; cfg.data.fsample = 200;
        cfg.data.time = mat2cell(0:1/cfg.data.fsample:(length(cfg.data.trial)*...
            length(cfg.data.time)/cfg.data.fsample)-(1/cfg.data.fsample), 1,length(cfg.data.trial)*length(cfg.data.time));
        cfg.data=rmfield(cfg.data,[{'sampleinfo'} {'trialinfo'}]);  % we don't need them inside the function to be able to filter continuous with no issues
        [data.so_measures] = lv_detect_phenomenon(cfg);
        % getting trough info.
        SO_trough_pos = zeros(1,length(cfg.data.time{1, 1}));
        for i=1:size(data.so_measures.measures,1), [~,id]=min(data.so_measures.measures{i,2});
            pos = data.so_measures.measures{i,1}(id);
            SO_trough_pos(pos) = 1;
        end
        SO_trough_mx = SO_trough_pos;
        SO_trough_mx = (reshape(SO_trough_mx,[sz(3) sz(1)]))';

%         good_trials = SO_trough_mx;
%         return; % to lock to trough of the SO
        %spindles
        cfg=[];
        cfg.phenomenon = 'spindle';
        cfg.fun_seq = [{'band'}  {'spec_thresholded'}];
        cfg.data = dat; cfg.data.fsample = 200;
        cfg.data.time = mat2cell(0:1/cfg.data.fsample:(length(cfg.data.trial)*...
            length(cfg.data.time)/cfg.data.fsample)-(1/cfg.data.fsample), 1,length(cfg.data.trial)*length(cfg.data.time));
        cfg.data=rmfield(cfg.data,[{'sampleinfo'} {'trialinfo'}]);  % we don't need them inside the function to be able to filter continuous with no issues
        [data.spindle_measures] = lv_detect_phenomenon(cfg);

        spindle_pos = zeros(1,length(cfg.data.time{1, 1}));
        for i=1:size(data.spindle_measures.measures,1)
            if isfield(cfg_hold,'lock_to_spindle_trough')
                [~,troughs] = min( data.spindle_measures.measures{i, 2});
                % troughs = islocalmin( data.spindle_measures.measures{i, 2});
                % putting one before the trough 50ms position so that it
                % becomes the center so we see the 50ms before it and the 50ms after
                cut_before = 25; % cut after is when you lock to event
                %                 spindle_pos(data.spindle_measures.measures{i,1}(troughs) - cfg.data.fsample*(cut_before/1000)) = 1;
                if round(mean(data.spindle_measures.measures{i,1})) - cfg.data.fsample*(cut_before/1000) >0 % to make sure wee don't go back before the first idx
%                     spindle_pos(round(mean(data.spindle_measures.measures{i,1})) - cfg.data.fsample*(cut_before/1000)) = 1; % spindle center
                    
                    spindle_pos( round(mean(data.spindle_measures.measures{i,1})) - cfg.data.fsample*(cut_before/1000) ...
                        : round(mean(data.spindle_measures.measures{i,1})) + cfg.data.fsample*(cut_before/1000) ) = 1; % spindle center

                end
            else
                spindle_pos(data.spindle_measures.measures{i,1}) = 1; % spindle_pos(data.spindle_measures.measures{i,1}(1)) = 1; for the beginning of the spindle only
            end
        end

        spindle_mx = spindle_pos; % take indivdual spindles and look at them they look good
        spindle_mx = (reshape(spindle_mx,[sz(3) sz(1)]))';
        if isfield(cfg_hold,'lock_to_spindle') % locking to spindles or spindles' beginning
            good_trials = spindle_mx; return;
        end
        good_trials = zeros(size(spindle_mx));
        % now we look at each spindle and see if there was a SO trough
        % before it in the same trial .. we look at 2.5 back in time from spindle beginning
        distance_time = 2.5; distance_samples = (distance_time*cfg.data.fsample);
        for i=1:size(spindle_mx,1)
            if ~isempty(find(spindle_mx(i,:)==1)) && ~isempty(find(SO_trough_mx(i,:)==1))
                y = find(SO_trough_mx(i,:)==1);
                spindle_SO_distance = cell2mat(arrayfun(@(x) (x-y), find(spindle_mx(i,:)==1) ,'Un',0)); % calculates the distance from every spindle beginning to every SO trough in the same trial

                if any(spindle_SO_distance <= distance_samples) && ...
                        any(spindle_SO_distance > 0) % because we look for spindles after trough
                    good_trials(i,:) =  SO_trough_mx(i,:); %spindle_mx(i,:); SO_trough_mx(i,:); % now good_trials has the spindle with SO before it in the same trial
                end
            end
        end
    case 'PAC' % phase amplitude coupling
        %% actual coupling
        % get SO and spindles (or any two phenom.) and then get the hilbert
        % analytical signals to know the complex info during the time of the events
        % then perform the actual coupling by taking the magnitude of the fast
        % phenom. and the phase of the slower and put them in euler's form and then
        % correct the result with surrogate data to have a modulation index which
        % is complex valued: its magnitude reflects the coupling strength and its
        % phase, the (slower phenom) phase for which there is the highest (faster phenom) amplitudes
        % in this example we use SO as the slower phenom and spindles as the faster
        % phenom.. if we want another phenom just consider the slower as so and the faster
        % as spindle so you change the band of detecting the phenom and continue
        % with the updated bands
        % you may construct a vector of records were reactivation happened and
        % check the coupling strength of that vs incorrect or control beaware of
        % the channels used for SO and spindles ..
        % preparing data into the continuous format
        if isfield(cfg,'keepidx'), keepidx = reshape( squeeze(cfg.keepidx(:,ismember(cfg.data.label,'fz'),:))' ,[],1);
            event_amp = keepidx; % that's posterior prob. or another amplitude other than spindle for which we do the PAC
            keepidx = find(~isnan(keepidx)); % keep the non nans indices
         end 

        dat = cfg.data; sz = size(dat.trial);
        dat.sampleinfo = (1:length(dat.time):length(dat.time)*sz(1))';
        dat.sampleinfo = [dat.sampleinfo dat.sampleinfo+length(dat.time)-1];
        cfg = []; cfg.channel = 'fz'; warning(['channel of phenomenon is set to: ' cfg.channel]);
        dat = ft_selectdata(cfg,dat);
        temp = squeeze(dat.trial); dat.trial=[];
        for i=1:size(temp,1), dat.trial{1,i}=temp(i,:); end % ch_time

        % formatting and SO
        cfg=[];
        cfg.phenomenon = 'SO';
        cfg.fun_seq = [{'band'} {'duration_zcrossings'} {'spec_thresholded'}];
        cfg.data = dat; cfg.data.fsample = 200;
        cfg.data.time = mat2cell(0:1/cfg.data.fsample:(length(cfg.data.trial)*...
            length(cfg.data.time)/cfg.data.fsample)-(1/cfg.data.fsample), 1,length(cfg.data.trial)*length(cfg.data.time));
        cfg.data=rmfield(cfg.data,[{'sampleinfo'} {'trialinfo'}]);  % we don't need them inside the function to be able to filter continuous with no issues
        data = cfg.data; % data is now formatted and can be used

        cfg.data.fs = 200;
        cfg.data = data;
        [so_measures] = lv_detect_phenomenon(cfg); % to get the timing of SO

        % those commented below should go with max_spindle .. because
        % cutting is done here but for the other cases it is done inside
        % lv_detect_phenomenon
%         if exist('keepidx','var'), temp = find(cell2mat(cellfun(@(x) ~isempty(intersect(x,keepidx)),so_measures.measures(:,1),'Un',0))); 
%             so_measures.measures = so_measures.measures(temp,:); so_measures.idx = so_measures.idx(temp,:);
%         end
%         so_measures.measures(:,1) = (cellfun(@(x) (intersect(x,keepidx))',so_measures.measures(:,1),'Un',0)); % for keeping exactly the indices which 
        % is needed in the correct pts analysis
        
        %spindles
        cfg=[];
        cfg.phenomenon = 'spindle';
        cfg.fun_seq = [{'band'}  {'spec_thresholded'}];
        cfg.data.fs = 200;
        cfg.data = data;
        [spindle_measures] = lv_detect_phenomenon(cfg); % to get the timing of spindles
        
        if exist('keepidx','var'), temp = find(cell2mat(cellfun(@(x) ~isempty(intersect(x,keepidx)),spindle_measures.measures(:,1),'Un',0))); 
            spindle_measures.measures = spindle_measures.measures(temp,:); spindle_measures.idx = spindle_measures.idx(temp,:);
        end
        % SO_hilbert
        cfg=[];
        cfg.phenomenon = 'Instant analytical';
        cfg.fun_seq = [{'band'}  {'spec_thresholded'}]; %0.16 1.25
        cfg.specs = [{'Instant analytical'} {'0.5 2'} {''} {'instant_phase_magnitude'} {'no'}]; % override the record with new values
        cfg.data=data;
        so_hilbert= lv_detect_phenomenon(cfg); % to get the phase of SO
        % Spindle_hilbert
        cfg=[];
        cfg.phenomenon = 'Instant analytical';
        cfg.fun_seq = [{'band'}  {'spec_thresholded'}]; % 12 16
        cfg.specs = [{'Instant analytical'} {'12 18'} {''} {'instant_phase_magnitude'} {'no'}]; % override the record with new values
        cfg.data=data;
        spindle_hilbert = lv_detect_phenomenon(cfg); % to get the magnitude of spindles

        % changing the amplitude of PAC to something other than spindle (e.g., certainty)
        spindle_measures.measures=[]; spindle_measures.measures{1,1} = keepidx;
        spindle_hilbert=[]; spindle_hilbert.measures{1, 3}=event_amp;

        % Coupling
        cfg=[];
        cfg.phenomenon = 'SO_spindles complexes';
        cfg.fun_seq = {'coupled_SOspindle'};
        cfg.data=data; % just to see data but with coupling we use the extracted measures
        cfg.data.so_measures = so_measures;
        cfg.data.spindle_measures = spindle_measures;
        cfg.data.so_hilbert = so_hilbert;
        cfg.data.spindle_hilbert= spindle_hilbert; % use the idx_limiter to choose a toi forr which PAC is calculated
        good_trials = lv_detect_phenomenon(cfg); % good_trials here got: coupling_strength and phase_of_coupling

    case 'fidelity'
        % returns a flag matrix of good_trials marking the beginning of the
        % most certain time pts according to thresholding criteria to use
        % them later in time locking to the most certain pt. or percentile
        % in case we want reactivation might be recurrent
        % the input is dval which comes from classification so they are the
        % posterior or z=w'*x values... and then we apply the thresholding
        % here and return the flags matrix ... either one most certain per
        % trial or percentile ...
        dval = cfg.data; good_trials=zeros(size(dval));% trnpts_trl_time
        if isfield(cfg,'percentile') 
            sz = size(dval);
            good_trials = nan(sz); c=0;
            % local maxima of every cerainty across time
            for i=1:sz(1)
                for j=1:sz(2)
                    [pks,locs] = findpeaks(squeeze(dval(i,j,:)), 'MinPeakDistance',0.1*200,...
                        'NPeaks',5,'SortStr','descend'); % 100ms*samplingRate
                    good_trials(i,j,locs) = pks; c=c+1;
                    lv_progress(c,sz(1)*sz(2), 'getting peaks');
                end
            end
            temp = good_trials(:);
            threshold=prctile(temp(~isnan(temp)),95)
            temp(temp>threshold)=1; temp(temp~=1)=0;
            good_trials = reshape(temp, sz);
        else % max certainty pt from each trial
            [~,id] = max(dval,[],3);
            for i=1:size(dval,1)
                for j=1:size(dval,2)
                    good_trials(i,j,id(i,j)) = 1;
                end
            end
        end




    otherwise
        fprintf(['\n Please select a method! \n']); good_trials = [];
end


end








% helping function
function [mu_rho] = predict_corr(trn, tst)
% takes train data and test data 2d(trl_feat)
% and returns the mean correlation between each of the test trl and all
% training trls

rho = corr(tst', trn', 'type','spearman'); % ex: 1st row the frst col from sleep corr with all trls from wake
mu_rho = mean(mean(rho, 2));

end




