function features_data = lv_feature_extractor(cfg)
% takes struct with trial (trl_ch_time), label and time and returns the same struct but
% now with features not normal EEG

if isfield(cfg,'data')
    data = cfg.data;
    if ~isfield(cfg,'method'), cfg.method = 'mean'; cfg.window=80; end
    if strcmp(cfg.method,'erp')==1, cfg.method='mean';
        if ~isfield(cfg,'fs') % uncomment for time domain smoothing
            sampling_rate = length(nearest(data.time,0):nearest(data.time,1)) - 1; warning(['sampling rate is set to: ' num2str(sampling_rate)]);
        else
            sampling_rate = cfg.fs; warning(['sampling rate is set to: ' num2str(sampling_rate)]);
        end
    end

    features_data = data;
end


switch cfg.method
    case {'phase','power'} % hilbert at every time pt .. be careful of edges
        warning('be careful of edges');
        hilbert_dat=[];
        cfg_preprocessing                 = [];
        cfg_preprocessing.bpfilter        = 'yes';
        %range = lv_tune_params('frequency band','4 8');
        cfg_preprocessing.bpfreq          = [9 13]; %str2double(split(range)');
        data_bp= ft_preprocessing(cfg_preprocessing, data);
        % matlab's hilbert
        for h=1:size(data_bp.trial,1) %trls
            for j=1:size(data_bp.trial,2) %ch
                hilbert_dat(h,j,:) = hilbert( squeeze(data_bp.trial(h,j,:)) ); % time should be in the first dimension.
            end
        end % power: is the squared magnitude of the complex vector at an instant in time. ... abs(hilbert_dat).^2 ... phase: angle(hilbert_dat)
        if strcmp(cfg.method,'power')==1
            features_data.trial = abs(hilbert_dat).^2;
        else
            phase_vals = angle(hilbert_dat);
            features_data.trial = [cos(phase_vals) sin(phase_vals)];
        end

    case 'spectrum'
        % the third dim should be time .. keeps the dims and get the
        % fourier coeffs of the signals
        dat = data.trial;
        nyquistVar = cfg.sampling_rate/2; N = size(dat,3);
        pts_hz = linspace(0,nyquistVar,(N/2)+1);
        coeff = fft(dat,[],3);
        features_data.trial = (abs(coeff(:,:,1:length(pts_hz)))*2)/N;  % spectrums

    case 'time_window' % in this one we don't take the mean of the time pts. we put them together in the FV (in channel dimension)
        wind_time = cfg.wind_time;
        fsample = data.fsample;
        timeax = data.time;

        % 20 samples = 100ms
        sample_time = 1*(1000)/fsample; %one sample in time
        temp=[]; temp_time=[];
        samples = (wind_time / sample_time)+1;  %41 200ms;

        samples = length(nearest(timeax,0):nearest(timeax,wind_time/1000));

        for i=1:samples
            temp = [ temp data.trial(:,:,i:end-(samples-i))];
            temp_time = [temp_time mean(timeax(i:i+samples-1))];
        end
        features_data.trial = temp;
        features_data.time = temp_time;  %timeax(ceil(samples/2):end-floor(samples/2));

    case {'log_var'} % when you think there is increased var/pw in specific band .. time is compressed .. because like ion motor imagery power
        % increase/decrease is not locked to stim it's not erp component it's self paced asynchronous event
        warning('be careful of edges');
        cfg_preprocessing                 = [];
        cfg_preprocessing.bpfilter        = 'yes';
        cfg_preprocessing.bpfreq          = [8 12];
        data_bp= ft_preprocessing(cfg_preprocessing, data);

        features_data= data;
        features_data.trial = log(var(data_bp.trial,[],3));

    case 'dwt_decomposition'
        % wavelet decomposition and in this case we are working with
        % signals and the input should be 1d and it has another version
        % wavedec2 which works on images so can be used on time freq. plots
        % or images.. downsamples the signal to half the points with every level
        % keeps some temporal and spectral features..
        cfg.levels = 3; % to keep the length of the signal and not increase it we have this no. of levels which could be changed 
        features_data= data;
        for i=1:size(data.trial,1)
            for j=1:size(data.trial,2)
                [c,L] = wavedec(squeeze(data.trial(i,j,:)),cfg.levels,cfg.wavelet);
                features_data.trial(i,j,:) = [appcoef(c,L,cfg.wavelet) detcoef(c,L,[1 2 3])]; % approximation of lvl3 then details of lvl1 then 2 then 3
            end
        end

    case {'dist_from_mean','dist_from_mean_cov'}
        mean_func = str2func(cfg.function);
        id1 = find(data.trialinfo==1); id2 = find(data.trialinfo==2);
        mu1 = mean_func(data.trial(id1,:,:),1);
        mu2 = mean_func(data.trial(id2,:,:),1);

        if strcmp(cfg.method,'dist_from_mean')==1
            y1=repmat(mu1,size(data.trial,1),1,1 ); y2=repmat(mu2,size(data.trial,1),1,1 );
            features_data= data;
            features_data.trial = cat(2, abs(data.trial-y1),abs(data.trial-y2));
        else % cov matrix for every trial and then distance from that
            for i=1:size(data.trial,1), mx(i,:) = reshape(cov((squeeze(data.trial(i,:,:) ))'), 1,[]); end % reshaping the cov matrix into vector
            features_data= data;
            features_data.trial = [mx - repmat( mean(mx(id1,:),1),size(data.trial,1),1 ) ...
                mx - repmat( mean(mx(id2,:),1),size(data.trial,1),1 ) ];
            % distance from mean covariance of first class concatenated
            % with the dist. for the second class .. prodeces
            % trialsx(ch_ch*2) because ch_ch is the cov put as vector and
            % then the two classes means distances is the 2

            % covert it to upper of the covariance because it's
            % symmetric!
        end
    case 'encoded_pts'
        % every point of every trial and channel and time will be encoded as
        % a 4d point with value channelx channely time ..
        sz=size(data.trial);
        times=1:sz(3); features_data.trial=[];
        if(~isequal(data.layout.label,data.label)), error('channel mismatch!'); end

        temp(1,:,:)=repmat(times,sz(2),1); times = repmat(temp,sz(1),1,1); clear temp; % time
        temp(1,:,:)=repmat(data.layout.pos(:,1),1,sz(3)); chx = repmat(temp,sz(1),1,1); clear temp; % channels
        temp(1,:,:)=repmat(data.layout.pos(:,2),1,sz(3)); chy = repmat(temp,sz(1),1,1); clear temp; % channels
        temp(:,1,:)=repmat(data.trialinfo(:,1),1,sz(3)); features_data.trialinfo=repmat(temp,1,sz(2),1); clear temp;
        features_data.trial = [data.trial(:) times(:) chx(:) chy(:)];
        features_data.trialinfo = features_data.trialinfo(:);

    case 'spatial_smoothing'
        % smoothing the data spatially using fieldtrip's neighbors function
        cfg_hold = cfg;
        cfg = []; cfg.method = 'triangulation'; cfg.senstype = 'EEG';
        cfg.layout = cfg_hold.layout; cfg.feedback = 'no';
        neighbours	= ft_prepare_neighbours(cfg);

        features_data = data;
        for i=1:length(neighbours)
            temp = neighbours(i); idx=ismember(cfg_hold.data.label, temp.label);
            idx_neighbors=ismember(cfg_hold.data.label, temp.neighblabel);
            features_data.trial(:,idx,:) = (0.8.*features_data.trial(:,idx,:) + 0.2.*mean( data.trial(:,idx_neighbors,:),2 ));
        end
    case 'crossfreq_PAC'
        % takes data and two frequency ranges and gets the PAC inside of
        % those ranges with pval
        warning('be careful of edges and give this function one channel in trl_ch_time');
        method = 1; % hilbert or fourier.. with hilbert the freq. are compressed but with fourier it is 2d plot
        if method==1
            hilbert_dat_low=[]; hilbert_dat_high=[];
            cfg_preprocessing                 = [];
            cfg_preprocessing.bpfilter        = 'yes';
            cfg_preprocessing.bpfreq          = cfg.freqLow;
            data_bp_low= ft_preprocessing(cfg_preprocessing, data);
            cfg_preprocessing.bpfreq          = cfg.freqHigh;
            data_bp_high= ft_preprocessing(cfg_preprocessing, data);

            cfg_temp=[]; cfg_temp.latency=cfg.toi_window; data_bp_low=ft_selectdata(cfg_temp,data_bp_low);
            data_bp_high=ft_selectdata(cfg_temp,data_bp_high);
            % matlab's hilbert
            if ~isequal(size(data_bp_low.trial),size(data_bp_high.trial)), error('sizes mismatch between data of low and high frequencies'); end
            for h=1:size(data_bp_low.trial,1) %trls
                for j=1:size(data_bp_low.trial,2) %ch
                    hilbert_dat_low(h,j,:) = hilbert( squeeze(data_bp_low.trial(h,j,:)) ); % time should be in the first dimension.
                    hilbert_dat_high(h,j,:) = hilbert( squeeze(data_bp_high.trial(h,j,:)) );
                end
            end % power: is the squared magnitude of the complex vector at an instant in time. ... abs(hilbert_dat).^2 ... phase: angle(hilbert_dat)

            amplitude = abs(hilbert_dat_high);
            phase = angle(hilbert_dat_low);
            amplitude = reshape(squeeze(amplitude)', 1,[]);    % all trials will be aggregated together in trial_time vector
            phase = reshape(squeeze(phase)', 1,[]);
        else
            freqLow = cfg.freqLow; freqHigh = cfg.freqHigh; toi_window = cfg.toi_window;
            cfg              = [];
            cfg.output       = 'fourier';
            cfg.channel      = 'all';
            cfg.method       = 'mtmconvol';
            cfg.taper        = 'hanning';
            cfg.toi          = data.time; % .time for max resolution .. window(1):0.1:window(2) the jumps just for visual smoothing
            cfg.pad          ='nextpow2'; % faster estimation
            
            frequencies = freqLow;
            cfg.foi          = linspace(frequencies(1),frequencies(end),2*(1+frequencies(end)-frequencies(1)));% frequencies(1):0.5:frequencies(end) default: 2 to 30 HZ
            % foi is arbitrary if big then it will repeat values but still correct
            cfg.t_ftimwin    =  3./cfg.foi; %not dynamic because we don't want to have nans
            TFdat_low = ft_freqanalysis(cfg, data);
            cfg_temp=[]; cfg_temp.latency=toi_window; TFdat_low=ft_selectdata(cfg_temp,TFdat_low);
            
            frequencies = freqHigh;
            cfg.foi          = linspace(frequencies(1),frequencies(end),2*(1+frequencies(end)-frequencies(1)));
            cfg.t_ftimwin    =3./cfg.foi; % 0.5*ones(length(cfg.foi),1); %
            TFdat_high = ft_freqanalysis(cfg, data);
            cfg_temp=[]; cfg_temp.latency=toi_window; TFdat_high=ft_selectdata(cfg_temp,TFdat_high);
            
            TFdat_high.fourierspctrm = reshape( permute( squeeze(TFdat_high.fourierspctrm),[2 3 1] ), length(TFdat_high.freq),[]);
            TFdat_low.fourierspctrm = reshape( permute( squeeze(TFdat_low.fourierspctrm),[2 3 1] ), length(TFdat_low.freq),[]);

            amplitude = abs(TFdat_high.fourierspctrm); % freq_time
            phase = angle(TFdat_low.fourierspctrm);
            result.hfreq = TFdat_high.freq;
            result.lfreq = TFdat_low.freq;
        end

        
        % phase/amplitude is a matrix with each row representing a freq. phase/amplitude
        option = 1; % 1 for Canolty's method and 2 for Tort's

        if option==1
            [result.coupling_strength, result.phase_of_coupling] = data2mvl(phase,amplitude); % this reshape is in case of one row
            nbin=[];
        elseif option==2
            nbin= 20;
            [result.coupling_strength, result.phase_of_coupling] = data2pac(phase,amplitude, nbin); % 20 is the number of bins (in tort's method they set it to 18 and ft to 20)
        end
        % Permutation to get the pvalue of coupling strength and know if it's significant
        % the permuted amplitude-phase assignment is done via cutting the amplitude in half at random pts. and re-run
        permutations = 1000; rng(1);
        sp = randi([2,size(phase,2)-1], permutations,1); % random splits for amplitude to create new amplitude-phase data
        parfor i=1:size(sp,1)
            amplitude_temp = [amplitude(:,sp(i,1):end) amplitude(:,1:sp(i,1)-1)];
            if option==1
                coupling_strength_temp(:,:,i) = data2mvl(phase,amplitude_temp); % this reshape is in case of one row
            elseif option==2
                coupling_strength_temp(:,:,i) = data2pac(phase,amplitude_temp, nbin); % 20 is the number of bins in tort's method
            end
        end
        result.permutations = coupling_strength_temp;
        
        % Pval
        for i=1:size(result.coupling_strength,1)
            for j=1:size(result.coupling_strength,2)
                result.pval(i,j) = mean(squeeze(coupling_strength_temp(i,j,:)) > result.coupling_strength(i,j));
            end
        end
        features_data = result;

    case 'topo_to_hilbert_curve'
        % converts the channels proximity from 2d to 1d vector and put it
        % back in the place of channels .. it's good because the number
        % will be high for channels because we got many pixels and the
        % spatial features are accurate and any jitter in them is captured in the curve in a smooth way
        lv_layout = data.layout; data.trial = single(data.trial);
        dim = 64;
        img = reshape(1:dim*dim, dim,dim);
        [hilb_idx  ,path] = hilbertCurve(img); % because we are feeding the indices in img we will get the indices
        temp =  (zeros(size(data.trial,1),dim*dim,size(data.trial,3)));
        for i=1:size(data.trial,1)
            lv_progress(i,size(data.trial,1),'progress: ');
            for j=1:size(data.trial,3)
                [topo, ~]=ft_plot_topo(lv_layout.pos(:,1), lv_layout.pos(:,2), squeeze(data.trial(i,:,j)),...
                    'mask', lv_layout.mask,'outline' ,lv_layout.outline);
                topo = imresize(topo,[dim dim]);% new dim 64*64 will give 4096 features
                temp(i,:,j) =topo(hilb_idx);
            end
        end
        features_data.trial = temp;
        figure
        imagesc(topo); set(gca,'YDir','normal'); hold on,
        plot(path(:,1),path(:,2))
    case 'hilbert_curve_to_topo'
        % converts a hilbert curve in the third dimension (like time) to a
        % 2d image to be fed to a classifier like CNN that takes advantage
        % of images and spatial proximity in 2d
        lv_layout = data.layout;
        % interpolating channels locations to fit in many subplots on the screen
        pos = [interp1([min(lv_layout.pos(:,1)) max(lv_layout.pos(:,1))],[0.1 0.9],lv_layout.pos(:,1)) ...
            interp1([min(lv_layout.pos(:,2)) max(lv_layout.pos(:,2))],[0.1 0.9],lv_layout.pos(:,2))];

        dim = 64;
        for i=1:size(data.trial,1)
            for j=1:size(data.trial,2)
                axes('Position',[pos(j,1) pos(j,2)  0.1 0.1]) % awel 2 el x w el y aken el figure nafsaha mn 0 l7ad 1 w da el makan elle feh el figure
                temp = imresize(squeeze(data.trial(i,j,:))',[1 dim*dim]);
                [hilb_img ,path] = hilbertCurve(temp); % from vector to image
                imagesc(hilb_img);
            end
            s=1; % to see the image
            figure(1)
            F = getframe(gcf);
            X = frame2im(F); %imagesc(X)
            features_data.trial{i,1}=X; % the image of every trial
        end
    case 'eeg_to_topos_video'
        % converts EEG to sequence of topos at every time pt. it's not an
        % actual video but instead of having the channels in the second
        % dimension we will have the image of the topo reshaped into a
        % vector and then we reshape back to 2d in a classifier like CNN
        % that uses the images in classification for 2d proximity
        lv_layout = data.layout; data.trial = single(data.trial);
        %         dim = 100; % was 10
        %         temp = nan(size(data.trial,1),dim*dim, size(data.trial,3));
        %         progressbar = ParforProgressbar(size(data.trial,1),'title', 'EEG to topo figures, progress');
        sz=size(data.trial,3);
        parfor i=1:size(data.trial,1)
            for j=1:sz
                %                                 [topo, ~]=ft_plot_topo(lv_layout.pos(:,1), lv_layout.pos(:,2), squeeze(data.trial(i,:,j)),...
                %                                     'mask', lv_layout.mask,'outline' ,lv_layout.outline);
                topo = topomaps(lv_layout,squeeze(data.trial(i,:,j)));
                topo = fillmissing(topo,'nearest'); topo(:,1)=[]; topo(:,end)=[];
                topo = imresize(topo,[100 100]);
                temp(i,:,j) = reshape(topo, 1,[]);% was reshape(imresize(topo,[dim dim]), 1,[]);
            end
            %             progressbar.increment();
        end
        %         delete(progressbar);
        close all;
        features_data.trial = temp;
    case 'CNN'
        % converts images to vectors using the activations of the last
        % layer so that it uses the proximity in the 2d image and converts
        % that to a meaningful vector .. that can be use with RNN or
        % another classifier the nice thing is that it should capture the
        % spatial features
        if cfg.transform==0
            data = data.trial; labels = cfg.data.trialinfo;
            % reshaping back from a vector in the second dimension to the 2d topo at the 2nd dim that it should be
            % the dimesions of topos images are fixed to be: 64*64
            data = reshape(data, size(data,1),sqrt(size(data,2)),sqrt(size(data,2)) );
            temp(1,:,:,:) = data;
            % numeric array of images should have this order of dimensions (hight width no.channels(color) no.images)
            data = permute(temp,[3 4 1 2]); clear temp;

            labels=categorical(labels);

            layers = [
                imageInputLayer([size(data,1) size(data,1) 1])

                convolution2dLayer(3,8,'Padding','same') % filter_size = 3x3 .. no.filters = 8
                batchNormalizationLayer
                reluLayer

                maxPooling2dLayer(2,'Stride',2) % downsampling 2x2 maximum inside this mask so it's downsampling ..
                % stride is the jumps so we skip a point we don't slide at every pt

                convolution2dLayer(3,16,'Padding','same')
                batchNormalizationLayer
                reluLayer

                maxPooling2dLayer(2,'Stride',2)

                convolution2dLayer(3,32,'Padding','same')
                batchNormalizationLayer
                reluLayer

                maxPooling2dLayer(10,'Stride',10)   % added this one to get the result with a more reduced dimensionality


                fullyConnectedLayer(2)
                softmaxLayer
                classificationLayer
                ];

            options = trainingOptions('adam', ...
                'InitialLearnRate',0.01, ...
                'MaxEpochs',100, ...
                'Verbose',false); % 'Plots','training-progress' add this for plotting training performance

            net = trainNetwork(data,labels',layers,options);
            features_data = net;
        else
            net = cfg.net;
            % RNN transform
            dt = cfg.data_to_transform.trial;
            dt = reshape(dt, size(dt,1),sqrt(size(dt,2)),sqrt(size(dt,2)) );
            temp(1,:,:,:) = dt;
            dt = permute(temp,[3 4 1 2]); clear temp;
            % activation output for the data_to_transform
            layer_name = net.Layers(13, 1).Name; % because 12 is the one before the fully connected and then we feed this to another classifier like RNN
            features_data.trial = (activations(net, dt, layer_name,'OutputAs','columns' ) )';
        end
    case 'ch_time_to_hilbert3d'
        % converts a 3d space like ch_time to a 1d hilbert curve that
        % considers the spatial proximity
        lv_layout = data.layout; data.trial = single(data.trial);
        dim = 64;
        [x,y,z] = hilbert3(4);

        hilb_idx = round(interp1(linspace(-0.5,0.5,dim), 1:dim, [x' y' z'])); % idx must be rounded
        % yes this indexing of -0.5 and 0.5 can be improved because it's
        % not exactly 0.5 but it won't make a difference because we will
        % have the same number of points and the same curve anyway it won't make a difference in the order of features on the curve
        % visual test
        figure,
        plot3(hilb_idx(:,1),hilb_idx(:,2),hilb_idx(:,3))
        dat = rand(dim,dim,dim);
        pos = [62,7,7];
        dat(pos(1),pos(2),pos(3)) = 1000; hold on, scatter3(pos(1),pos(2),pos(3),'MarkerFaceColor','r') % displaying the position and the actual value in it is 1000
        for id=1:size(hilb_idx,1), temp(1,id) = dat(hilb_idx(id,1),hilb_idx(id,2),hilb_idx(id,3)); end
        figure, plot(temp) % to see where is 1000 in 1d

        temp = nan(size(data.trial,1),dim*dim);
        for i=1:size(data.trial,1)
            time_vec_query = linspace(min(data.time),max(data.time),dim);
            dat2d = (interp1(data.time, squeeze(data.trial(i,:,:))', time_vec_query))'; % ch_time
            lv_progress(i,size(data.trial,1),'progress: '); topo_img=[];
            for j=1:size(dat2d,2)
                [topo, ~]=ft_plot_topo(lv_layout.pos(:,1), lv_layout.pos(:,2), dat2d(:,j),...
                    'mask', lv_layout.mask,'outline' ,lv_layout.outline);
                topo_img = cat(3, topo_img, imresize(topo,[dim dim]));% 64*64 images of topos put after each other in time
            end
            for id=1:size(hilb_idx,1), temp(i,id) = topo_img(hilb_idx(id,1),hilb_idx(id,2),hilb_idx(id,3)); end
        end
    case 'ch_time_to_hilbert3d2'
        % the difference is that this one puts nans in the places of no
        % no channels instead of interpolating so it doesn't inflate the
        % data and should be faster .. so it will order the data according
        % to hilbert without changing it.. we fill the gaps with nans in
        % time and channels
        features_data.trial=[];
        lv_layout = data.layout; data.trial = single(data.trial);
        dim = pow2(ceil(log2(size(data.trial,3)))); % oversample time to the nearest ^2
        [x,y,z] = hilbert3(log2(dim)); % log2 is something i experimented with until the curve didn't miss any data points and that was log2(dim)

        hilb_idx = round(interp1(linspace(-0.5,0.5,dim), 1:dim, [x' y' z'])); % idx must be rounded

        for i=1:size(data.trial,1)
            data2d = nan(size(data.trial,2),dim);
            timeidx_highdim = round(interp1(linspace(min(data.time),max(data.time),dim),1:dim,data.time));
            data2d(:,timeidx_highdim) = squeeze(data.trial(i,:,:));

            lv_progress(i,size(data.trial,1),'progress: ');

            % interpolating to know the positions of electrodes in the higher dim figure
            xd = round(interp1( linspace(min(lv_layout.pos(:,1)),max(lv_layout.pos(:,1)), dim),1:dim, lv_layout.pos(:,1) ));
            yd = round(interp1( linspace(min(lv_layout.pos(:,2)),max(lv_layout.pos(:,2)), dim),1:dim, lv_layout.pos(:,2) ));

            topos = nan(dim,dim,size(data2d,2));

            for k=1:length(xd), topos(yd(k),xd(k),:)=data2d(k,:); end

            templin_topo= topos(:);
            temp = templin_topo( sub2ind([dim dim dim],hilb_idx(:,1),hilb_idx(:,2),hilb_idx(:,3)) );

            features_data.trial(i,:) = temp(~isnan(temp))';
        end
    case 'euclidean alignment'
        % EA takes data and aligns them using the inverse of the covariance matrix to give the precision matrix
        % which measures how close the data is to the mean how clustered they are and when that's
        % projected to the data it aligns it and then you can do across sbj classification and extract features from the aligned signals ..
        % following the paper: Transfer Learning for Brain�Computer Interfaces:A Euclidean Space Data Alignment Approach
        % equation 11 .. great for leave one sbj out classification
        if length(size(data.trial))==3, for i=1:size(data.trial,1), data_temp{1,i} = squeeze(data.trial(i,:,:)); end, data.trial=data_temp; end  % to 2d because easier with multiplication
        covs = pooled_cov(data.trial); inv_covs = chol(inv(covs)); % chol is Cholesky factorization which is the sqrt for the matrix to make inv(covs) = Result'*Result sqrt(x)*sqrt(x) = x

        proj_data = cellfun( @(x) ((inv_covs )*x), data.trial,'Un',0);

        % verify cov_projected: should be identity matrix
        cov_projected = pooled_cov(proj_data);

        % to 3d
        temp = cell2mat(proj_data);
        temp = reshape(temp,size(proj_data{1,1},1),size(proj_data{1,1},2),[]);
        features_data.trial = permute(temp, [3 1 2]);
    case 'euclidean alignment_trials' % for classification domain align between trn and tst not trials align
        % this is a manual edit but the paper is on trials of EEG as in the previous switch case
        source = cfg.trn.trial; % trials_features
        target = cfg.tst.trial;

        source_cov = cov(source); inv_cov = chol(inv(source_cov));
        proj_data = inv_cov*source'; features_data.source=proj_data';
        % eye matrix
        cov_projected = cov(features_data.source);

        target_cov = cov(target); inv_cov = chol(inv(target_cov));
        proj_data = inv_cov*target'; features_data.target=proj_data';

    case 'CORAL'
        % takes source and target and gives updated source, source is trn set data and aligns them using the source and target matrices
        % and takes the inverse of the covariance matrix of source and
        % whitenes it then we transform the other matrix with this whitened
        % cov it's a bit similar to EA but aligns two signals based on cov
        % matrices so it actually changes the data not just the alignment
        % we will use it here to align all signals inside data.trial like a
        % crosscorrelation alignment inside lv_align
        % data is assumed to be zscored .. this align is for the
        % rotation/covariance alignment
        source = cfg.trn; % trials_features
        target = cfg.tst;

        Cs = cov(source.trial) + eye(size(source.trial,2)); % identity with the no. features to make the cov full rank ... test:diag(var(target.trial,[],1)) when centered is identity
        Ct = cov(target.trial) + eye(size(target.trial,2));

        Ds = source.trial * chol(inv(Cs)); % whitening source ..

        features_data.source = Ds * chol(Ct); % re-coloring with target covariance to bring the source trials to the same feature space as the target
        % so that they become aligned in the feature space
        features_data.target = target.trial;
    case 'CORAL_manual'
        % same as coral but here we do some data visual inspection and by
        % looking at the projection I didn't do the eye summation and
        % recentered the data after projection.. because we need the source
        % projected in that direction directly...
        % data is assumed to be zscored .. this align is for the
        % rotation/covariance alignment
        source = cfg.trn; % trials_features
        target = cfg.tst;

        Cs = cov(source.trial); % identity with the no. features to make the cov full rank ... test:diag(var(target.trial,[],1)) when centered is identity
        Ct = cov(target.trial);

        Ds = source.trial * chol(inv(Cs)); % whitening source

        features_data.source = zscore(Ds * chol(Ct), [],1); % re-coloring with target covariance to bring the source trials to the same feature space as the target
        % so that they become aligned in the feature space
        features_data.target = target.trial;
    case 'subspace_alignment'
        % aligns the features of samples from two datasets the source and
        % the target by getting the pricipal component of the source and
        % the target and then transforming the PC of the source to fit the
        % target and by aligning the PCs and projecting the data to that
        % they align the data of source and target to a common space, from
        % Subsapce Alignment paper in methods papers
        % data is assumed to be zscored .. this align is for the
        % rotation/covariance alignment
        source = cfg.trn.trial; % trials_features
        target = cfg.tst.trial;

        [eigVects_source,~,eigVals] = pca(source); [~,id]=sort(eigVals,'descend');
        Xs = eigVects_source(:,id)' * source'; % eigVectors_trials
        Xs = Xs';
        %plot(Xs(:,1),Xs(:,2),'+')

        [eigVects,~,eigVals] = pca(target); [~,id]=sort(eigVals,'descend');
        Xt = eigVects(:,id)' * target'; Xt = Xt';

        M = eigVects_source(:,id)' * eigVects(:,id); % identity when the data of source and target are the same
        Xa = Xs*M; %trials_Vectors

        features_data.source = Xa; features_data.target = Xt;

    case 'subspace_alignment_manual'
        % testing aligning PCs using transformation matrix that maps
        % eigVectors of source to those of target Vt=R*Vs , so R = Vt*inv(Vs)and then getting this
        % transformation and applying it to data points
        % data is assumed to be zscored .. this align is for the
        % rotation/covariance alignment
        % in fact the sourcce is the testing set and trn is the target
        % becausse we rotate the tst to fit the trn
        source = cfg.trn.trial; % trials_features
        target = cfg.tst.trial;

        [eigVects_source] = pca(source); % the first has the highest variance pca already sorts based on the eigVals
        [eigVects_target] = pca(target);


        % manual pts of the PCs
        %         axis equal
        %         hold on, plot([0 eigVects_source(1,1)],[0 eigVects_source(2,1)]);
        %         plot([0 eigVects_source(1,2)],[0 eigVects_source(2,2)]);
        rot = eigVects_target' * (eigVects_source); % equals: rot = eigVects_target' * pinv(eigVects_source') but because that's eigVect so the inv is the transpose so we leave it
        proj = rot'*source';
        %         figure, plot(proj(1,1:100),proj(2,1:100),'+')
        %         hold on, plot(proj(1,101:200),proj(2,101:200),'o')
        features_data.source = proj'; features_data.target = target;
        features_data.source = zscore(features_data.source,[],1);

        %         eigVects_source = eigVects_source(:,1:2); eigVects_target = eigVects_target(:,1:2);
        %         features_data.source = zscore((eigVects_source' * source')', [],1);
        %         features_data.target = zscore((eigVects_target' * target')', [],1);

        % % %         % PCA
        % % %         features_data.source = (eigVects_target(:,1)' * source')';
        % % %         features_data.target = (eigVects_target(:,1)' * target')';

        % plot scattered points of original and projected
        %         subplot(121),plot(source(:,1),source(:,2),'+')
        %         hold on, plot(target(:,1),target(:,2),'o')
        %         subplot(122),plot(features_data.source(:,1),features_data.source(:,2),'+')
        %         hold on, plot(features_data.target(:,1),features_data.target(:,2),'o')
    case 'subspace_alignment_manual2'
        % each dataset will get its PCs and project on them
        source = cfg.trn.trial; % trials_features
        target = cfg.tst.trial;

        [eigVects_source] = pca(source); % the first has the highest variance pca already sorts based on the eigVals
        [eigVects_target] = pca(target);

        target_comp = 4;
        eigVects_target = eigVects_target(:,1:target_comp);
        eigVects_source = eigVects_source(:,1:target_comp);

        target = (eigVects_target' * target')';
        source = (eigVects_source' * source')';

        features_data.source = source; features_data.target = target;
    case 'ssa'
        % this one is on EEG .. it finds the stationary and non-stationary sources of the signals ..
        % the stationary source is the one that maximises the dependence
        % between the signals and non stationary is noise.. guess the
        % signals should be carrying the same effect so within class..
        [data,data2]=deal([]);
        for i=1:size(cfg.trn.trial,1), data{1,i} = reshape(cfg.trn.trial(i,:,:),size(cfg.trn.trial,2),size(cfg.trn.trial,3)); end
        for i=1:size(cfg.tst.trial,1), data2{1,i} = reshape(cfg.tst.trial(i,:,:),size(cfg.tst.trial,2),size(cfg.tst.trial,3)); end
        sources=1; % retrieve one source
        sz = size(cfg.trn.trial); sz2 = size(cfg.tst.trial);
        [est_Ps, est_Pn, est_As, est_An, ssa_results] = ssa(cell2mat(data)...
            , sources,'reps', 50,  'equal_epochs', length(data), 'random_seed', 12345);

        src_trn = est_Ps * cell2mat(data);
        src_tst = est_Ps * cell2mat(data2);

        src = reshape(src_trn,[sources sz(3) sz(1)] );
        features_data.source_trn = permute(src,[3 1 2]);
        %ssa_results.n_src non stationary source (noise)
        src = reshape(src_tst,[sources sz2(3) sz2(1)] );
        features_data.source_tst = permute(src,[3 1 2]);
    case 'PCA_correlation'
        % correlating PCs of the source and target and leaving the PCs that
        % will have significant positive correlation meaning that they are
        % similar in the feature space and should transform each dataset in
        % its own scaling but in a common way
        source = cfg.trn.trial; % trials_features
        target = cfg.tst.trial;
        if cfg.pool==0 % pooling data together?
            %             rng(1), [eigVects_source,~,eigVals_source] = pca(source); % the first has the highest variance pca already sorts based on the eigVals
            rng(1), [eigVects_target,~,eigVals_target] = pca(target);
            %             [rho,pval] = corr(eigVects_source, eigVects_target);% src x target
            %             temp = (rho>0 & pval<0.01);
            %         source_comp = find(sum(temp,2)>0);
            %             target_comp = find(sum(temp,1)>0); %actual trn set

            % fixing the components
            %             target_comp = 4;
            % 80% of explained variance
            eigVals_target = eigVals_target./sum(eigVals_target);
            target_comp = find(cumsum(eigVals_target) > 0.9);
            target_comp = target_comp(1);
            % Kaiser-Guttman criterion
            %             target_comp = sum(eigVals_target>1);

            features_data.source = (eigVects_target(:,1:target_comp)' * source')';
            features_data.target = (eigVects_target(:,1:target_comp)' * target')';
        elseif cfg.pool==0.5 % dataset dependent transformation
            % max correlation and transform each dataset with its PCs
            rng(1), [eigVects_source] = pca(source); % the first has the highest variance pca already sorts based on the eigVals
            rng(1), [eigVects_target] = pca(target);
            [rho,pval] = corr(eigVects_source, eigVects_target);% src x target
            temp = (rho>0 & pval<0.01);
            masked_rho = double(temp).*rho;
            [val,source_id] = max(masked_rho,[],2); % be careful because the order matters because
            % we delete the col first then row and if we delete the row first it might be different
            for i=1:size(masked_rho,1), masked_rho(i, masked_rho(i,:)~=val(i))=0; end
            [val,target_id] = max(masked_rho,[],1);
            for i=1:size(masked_rho,2), masked_rho(masked_rho(:,i)~=val(i), i)=0; end
            [source_comp, target_comp] = find(masked_rho>0);

            features_data.source = (eigVects_source(:,source_comp)' * source')';
            features_data.target = (eigVects_target(:,target_comp)' * target')';
        elseif cfg.pool==1
            no_PCs = 4; %warning('using 4 PCs..');
            rng(1), [eigVects,~,eigVals_target] = pca([source;target]);
            % 80% of explained variance
            eigVals_target = eigVals_target./sum(eigVals_target);
            target_comp = find(cumsum(eigVals_target) >= 0.8);
            target_comp = target_comp(1);

            features_data.source = (eigVects(:,1:target_comp)' * source')';
            features_data.target = (eigVects(:,1:target_comp)' * target')';
        elseif cfg.pool==2 % kernel PCA for handling nonlinearities
            kpca = KernelPca(target, 'gaussian', 'gamma', 2.5, 'AutoScale', true);
            target_comp = 4;
            features_data.source = project(kpca, source, target_comp);
            features_data.target = project(kpca, target, target_comp);
        end
    case 'variance_featSelection'
        % removes features with lower than threshold% of variance because
        % features with low variance most likely don't contain useful info
        % will only work if data isn't zscored.. so we zscore after doing
        % this selection
        source = cfg.trn.trial; % trials_features
        target = cfg.tst.trial;

        var_target = var(target,[],1);
        [var_target, id] = sort(var_target,'descend');
        % 80% of explained variance
        var_target = var_target./sum(var_target);
        target_features = find(cumsum(var_target) >= 0.8);

        features_to_keep = id(1:target_features(1));
        features_data.source = zscore(source(:,features_to_keep),[],1);
        features_data.target = zscore(target(:,features_to_keep),[],1);

    case 'random_forest_featSelection'
        % feature selection using feature importance based on trees
        target = cfg.tst.trial; labels = cfg.labels;
        source = cfg.trn.trial;
        % ensemble of random forest
        t = templateTree('NumVariablesToSample','all',...
            'PredictorSelection','interaction-curvature','Surrogate','on');
        rng(1); % For reproducibility
        Mdl = fitcensemble(target,labels,'Method','Bag','NumLearningCycles',50, ...
            'Learners',t);
        impOOB = oobPermutedPredictorImportance(Mdl);
        [~,id] = sort(impOOB,'descend');
        features_data.source = source(:,id(1:4));
        features_data.target = target(:,id(1:4));
    case 'MI_featSelection'
        % Mutual Information based feature selection
        source = cfg.trn.trial;
        target = cfg.tst.trial; % actual trn set

        numToSelect=size(target,2); labels = cfg.labels;
        [selectedFeatures,scorespre,classMI]= MIM(numToSelect,target,labels);
        thresholdDistance = (scorespre(1)-scorespre(end)) / 10; % 10% of highest mutual info values
        sub = scorespre>scorespre(1)-thresholdDistance;
        features_data.target=target(:,selectedFeatures(sub));
        features_data.source=source(:,selectedFeatures(sub));
    case 'sequential_featSelction'
        % sequential_featSelction with SVM
        source = cfg.trn.trial;
        target = cfg.tst.trial; % actual trn set
        labels = cfg.labels; rng(1);
        c = cvpartition(labels,'k',5);
        opts = statset('UseParallel',true);
        fun = @(XT,yT,Xt,yt)loss(fitcdiscr(XT,yT),Xt,yt); % fitcdiscr:lda fitcecoc: multiclass SVM

        [fs,history] = sequentialfs(fun,target,labels,'cv',c,'options',opts);
        features_data.target=target(:,fs);
        features_data.source=source(:,fs);
    case 'short_time_cov'
        % takes the 3d data and slides on time with a window to get the
        % covariance for each ch_shorttime and then uses the upper triangle
        % of that cov as features .. good because it gets the change of
        % channels together by doing correlation that's what the cov does
        cfg2=[];
        cfg2.shift = cfg.shift; % shifts in samples between windows (minimum is one sample shift because then we will get the full resolution)
        cfg2.windlen = cfg.windlen; % in samples
        cfg2.data=data;
        result = lv_slider(cfg2); % samples_windows
        result.higher_data = single(result.higher_data);
        features_data.time = result.time;
        size(result.higher_data)
        upper=0; % to choose the upper triangle or the full cov

        if upper==1
            features_data.trial = single(nan(size(result.higher_data,1),(size(result.higher_data,2)*size(result.higher_data,2))-sum((1:size(result.higher_data,2))-1), ... % the sum is for the upper triangular
                size(result.higher_data,4) ));
            for i=1:size(result.higher_data,1)
                for w=1:size(result.higher_data,4)
                    temp = triu(cov( squeeze(result.higher_data(i,:,:,w))' )); % the upper trangular of the cov matrix
                    features_data.trial(i,:,w) = temp(temp~=0); % trial_features_window as if it's trial_ch_time and ch are features
                end
                lv_progress(i,size(result.higher_data,1),'cov''s upper is being put instead of channels: ');
            end
        else
            trial=[]; features_data.trial=[];
            sz = size(result.higher_data); t = sz(4);
            parfor i=1:sz(1)
                for w=1:t
                    temp = cov( squeeze(result.higher_data(i,:,:,w))' ); % full cov
                    trial(i,:,w) = temp(:);
                end
                %lv_progress(i,size(result.higher_data,1),'cov is being put instead of channels: ');
            end
            features_data.trial = trial;
        end
    case 'mdim_kernel'
        % takes 3d data and does smoothing but in 3d so it goes to
        % topo_time and then smoothes that and then gives the smoothed
        % output
        data.trial = single(data.trial);
        features_data.trial=[];
        lv_layout = data.layout; data.trial = single(data.trial);
        dim = 64; % the dimesion for the topo
        stride = 20;
        features_data.trial=[]; features_data.time=[];
        for i=1:size(data.trial,1)
            data2d = squeeze(data.trial(i,:,:));

            lv_progress(i,size(data.trial,1),'progress: ');

            % interpolating to know the positions of electrodes in the higher dim figure
            xd = round(interp1( linspace(min(lv_layout.pos(:,1)),max(lv_layout.pos(:,1)), dim),1:dim, lv_layout.pos(:,1) ));
            yd = round(interp1( linspace(min(lv_layout.pos(:,2)),max(lv_layout.pos(:,2)), dim),1:dim, lv_layout.pos(:,2) ));

            topos = single(nan(dim,dim,size(data2d,2)));

            for k=1:length(xd), topos(yd(k),xd(k),:)= data2d(k,:); end

            % imagesc(topos(:,:,10))  looking at this we can see that the
            % kernel should be as big as 20 in lengths o.w. the points won't see each other

            %             kernel = (1/(stride^3))*ones(stride,stride,stride);
            %
            %             convolved = convn(topos,kernel,'valid');
            %             for t=1:size(convolved,3), temp = convolved(:,:,t);
            %                 features_data.trial(i,:,t) = temp(:);
            %             end
            shift=1;  c1=0;
            winlen=stride;
            %
            %             if i==1 % getting the indices
            %                 idx = 1:numel(topos);
            %                 idx = reshape(idx, size(topos));
            %                 for d1=1:shift:size(topos,1)-(winlen-1)
            %                     c1=c1+1; c2=0; c3=0;
            %                     for d2=1:shift:size(topos,2)-(winlen-1)
            %                         c2=c2+1; c3=0;
            %                         for d3=1:shift:size(topos,3)-(winlen-1)
            %                             c3=c3+1; % perform the operation on the cubes here
            %                             cubes_cells{c1,c2,c3} = (idx(d1:d1+winlen-1, d2:d2+winlen-1, d3:d3+winlen-1));
            %                         end
            %                     end
            %                 end
            %             end
            %
            %             temp = topos(:); cubes_cells2 = cubes_cells(:);
            %             for j=1:numel(cubes_cells)
            %                 cubes_data{j} = temp(cubes_cells(j));
            %             end
            for d1=1:shift:size(topos,1)-(winlen-1)
                c1=c1+1; c2=0; c3=0;
                for d2=1:shift:size(topos,2)-(winlen-1)
                    c2=c2+1; c3=0;
                    for d3=1:shift:size(topos,3)-(winlen-1)
                        c3=c3+1; % perform the operation on the cubes here
                        temp= (topos(d1:d1+winlen-1, d2:d2+winlen-1, d3:d3+winlen-1));
                        cubes(c1,c2,c3) = nanmean(temp(:));
                        if i==1
                            times = repmat(data.time,size(topos,1),1);
                            times = permute(repmat(times,1,1,size(topos,2)), [1 3 2]);
                            tempt= (times(d1:d1+winlen-1, d2:d2+winlen-1, d3:d3+winlen-1));
                            features_data.time(c1,c2,c3) = nanmean(tempt(:));
                        end
                    end
                end
            end
            features_data.trial(i,:,:) = reshape(cubes,size(cubes,1)*size(cubes,2),[]);
            %             res = cell2mat(cellfun(@(x) (nanmean(nanmean(x,1),2)), cubes_data,'Un',0));
            %             features_data.trial(i,:,:) = reshape(res,size(res,1)*size(res,2),[]);
        end
        features_data.time = squeeze(mean(mean(features_data.time,1),2));

    case 'short_correlations' % short time correlations (difference to keep the direction) are correlations with a moving window
        % between lateral channels to get a correlation value across time
        % windows and classify them .. it takes the channels of two areas
        % and gets the correlations
        ch_left = cfg.ch_left; ch_right = cfg.ch_right;
        data = lv_slider(cfg);
        data.trial = data.higher_data;

        data_left = data.trial(:,cell2mat(cellfun(@(x) find(ismember(cfg.data.label,x)), ch_left,'Un',0)),:,:);
        data_right =  data.trial(:,cell2mat(cellfun(@(x) find(ismember(cfg.data.label,x)), ch_right,'Un',0)),:,:);
        features_data.trial=[];
        %         for i=1:size(data.trial,1), lv_progress(i,size(data.trial,1),'working on correlations: ')
        %             for ch=1:size(data_left,2)
        %                 rho = corr(squeeze(data_left(i,ch,:,:)),squeeze(data_right(i,ch,:,:)),'type','Spearman');
        %                 features_data.trial(i,ch,:) = diag(rho)';
        %             end
        %         end
        features_data.trial = squeeze(mean((data_left - data_right), 3));
        features_data.time = data.time;

    case 'multitapers_smoothing'
        sampling_rate = length(nearest(data.time,0):nearest(data.time,1)) - 1; warning(['sampling rate is set to: ' num2str(sampling_rate)]);
        window_samples = (sampling_rate*cfg.window)/1000;
        ws = floor(window_samples/2);
        time_idx = ws+1 : size(data.trial,3)-ws; % from the first possible window in time
        tapers = (dpss((2*ws)+1,3))'; 

        for tt=1:length(time_idx)
            temp = data.trial(:,:,time_idx(tt)-ws : time_idx(tt)+ws);
            for i=1:size(tapers,1)
                ss(1,:,:) = repmat(tapers(i,:),size(temp,2),1);
                tapered(:,:,i) = mean(bsxfun(@times, temp, ss) , 3); %trlxchxtapers
            end
            features_data.trial(:,:,time_idx(tt)) = mean(tapered ,3);
        end
        features_data.trial = features_data.trial(:,:,time_idx);
        features_data.time = features_data.time(1,time_idx);
 

    otherwise
        window_samples = (sampling_rate*cfg.window)/1000;
        ws = floor(window_samples/2);
        time_idx = ws+1 : size(data.trial,3)-ws; % from the first possible window in time
        for tt=1:length(time_idx)
            % any function that acts on a dimension as SECOND argument can be used: mean, median.. etc,
            features_data.trial(:,:,time_idx(tt)) = eval([cfg.method '( data.trial(:,:,time_idx(tt)-ws : time_idx(tt)+ws) ,3)']);
        end
        features_data.trial = features_data.trial(:,:,time_idx);
        features_data.time = features_data.time(1,time_idx);

end



end


% helping functions
function [x,y,z] = hilbert3(n)
% Hilbert 3D curve.
%
% function [x,y,z] = hilbert3(n) gives the vector coordinates of points
% in n-th order Hilbert curve of area 1.
%
% Example: plot the 3-rd order curve
%
% [x,y,z] = hilbert3(3); plot3(x,y,z)
%   Copyright (c) by Ivan Martynov
%   Inspired by function [x,y] = hilbert(n) made by Federico Forte
%   Date: September 17, 2009
if nargin ~= 1
    n = 2;
end

if n <= 0
    x = 0;
    y = 0;
    z = 0;
else
    [xo,yo,zo] = hilbert3(n-1);
    x = .5*[.5+zo .5+yo -.5+yo -.5-xo -.5-xo -.5-yo .5-yo .5+zo];
    y = .5*[.5+xo .5+zo .5+zo .5+yo -.5+yo -.5-zo -.5-zo -.5-xo];
    z = .5*[.5+yo -.5+xo -.5+xo .5-zo .5-zo -.5+xo -.5+xo .5-yo];
end
end


function cov_data_pooled=pooled_cov(data)
% gets the pooled covariance of different trials (trials in cells) or 3d
if length(size(data))==3, for i=1:size(data,1), data_temp{1,i} = squeeze(data(i,:,:)); end, data=data_temp; end % if 3d make it 2d
data_centered = cellfun(@(x) (x-repmat(mean(x,2),1,size(x,2)) ), data,'Un',0);
cov_data = cellfun(@(x) ((x*x')./size(x,2)), data_centered,'Un',0);
cov_data_pooled=zeros(size(cov_data{1,1}));
for i=1:length(data), cov_data_pooled = cov_data_pooled + cov_data{1,i}; end

cov_data_pooled = cov_data_pooled./length(data);

end



function topodata = topomaps(lv_layout,dat)
% from Mike's https://www.youtube.com/watch?v=VyAAW_fWo6M&t=1232s&ab_channel=MikeXCohen
% takes one vector and does interpoaltion with respect to the channels' locations
% topographical maps
% mikeXcohen@gmail.com

% get cartesian coordinates
elocsX = lv_layout.pos(:,1); elocsY = lv_layout.pos(:,2);

% define XY points for interpolation
interp_detail = 100;
interpX = linspace(min(elocsX),max(elocsX),interp_detail);
interpY = linspace(min(elocsY),max(elocsY),interp_detail);

% meshgrid is a function that creates 2D grid locations based on 1D inputs
[gridX,gridY] = meshgrid(interpX,interpY);

% now interpolate the data on a 2D grid
interpFunction = TriScatteredInterp(elocsX,elocsY,double(dat'), 'natural');
topodata = interpFunction(gridX,gridY);
topodata = single(topodata);
% figure, imagesc(topodata)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION from fieldtrip's crossfrequency for PAC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pacdata,lv_phase_angle] = data2mvl(LFsigtemp,HFsigtemp) % Canolty 2006
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

pacdata = abs(mvldata); % the magnitude of the modulation index is the coupling strength
lv_phase_angle = angle(mvldata); % angle is phase
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
for i=1:size(pacdata,1) %freqlow
    for j=1:size(pacdata,2) %freqhigh
        bins_temp = bins;
        w = squeeze(pacdata(i,j,:))./nansum(squeeze(pacdata(i,j,:)));
        bins_temp(isnan(w))=[]; w(isnan(w))=[];
        lv_phase_angle(i,j) = circ_mean(bins_temp, w);
    end
end

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