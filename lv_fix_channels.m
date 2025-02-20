function data_interp = lv_fix_channels(data)
choice = 0;

if choice==1
    % takes data struct and at every time point and every trial looks at the distribution among
    % the channels and interpolates the bad time points specified as outliers
    % the issue with this is that the variance may be high and still we get
    % the within trial stats. showing no outlier because all channels are
    % high together ..
    data.trial = single(data.trial);
    sz = size(data.trial);
    data.trial = permute(data.trial,[2 1 3]);
    data.trial = reshape(data.trial,sz(2),[])'; % trl_timexch to loop on trl_time only
    temp_time = data.time;
    data.time=[0 1]; lv_layout = data.lv_layout;
    goodFlags_matrix=nan(size(data.trial));
    parfor i=1:size(data.trial,1)
        all_measure = data.trial(i,:);
        iqr_val = iqr(all_measure); q1=prctile(all_measure,25); q3=prctile(all_measure,75);
        goodFlags_matrix(i,:) = double(all_measure>q1-(1.5*iqr_val) & all_measure<q3+(1.5*iqr_val)); % either is a bad channel
    end

    %% then we loop and put time pts in rows after trials and loop on that and then reshape ... so the time
    % dim will be one and the trials will be trial_time and then we loop then
    % reshape the interpolated result ...
    % Find neighbours
    cfg              = [];
    cfg.method       = 'triangulation'; % the way it's going to connect channels indicating neighborhood
    cfg.senstype     = 'EEG';
    cfg.layout       = lv_layout; cfg.feedback     = 'no';
    neighbours	= ft_prepare_neighbours(cfg);

    data.trial(:,:,2) = data.trial;
    data=rmfield(data,'trialinfo'); data=rmfield(data,'cfg'); data=rmfield(data,'lv_layout');

    tic
    data_interp = data;
    [unique_trls,~,ic] = unique(goodFlags_matrix,'rows');
    inter_co = 0; % interpolation count to make sure it was done correctly
    for i=1:size(unique_trls,1)
        idx = find(ic == i);

        cfg                  = [];
        cfg.badchannel       = data.label( ~goodFlags_matrix(idx(1),:) );
        if isempty(cfg.badchannel), data_interp.trial(idx,:,:) = data.trial(idx,:,:); continue; end % good trial
        cfg.neighbourdist    = 4;
        cfg.neighbours     	 = neighbours;
        cfg.layout           = lv_layout;
        cfg.method           = 'spline';
        cfg.feedback         = 'no';
        cfg.trackcallinfo	 = 'no';
        cfg.trials           = idx;

        interpAvg = ft_channelrepair(cfg ,data);
        if isfield(interpAvg,'avg'), interpAvg.avg = interpAvg.avg(:,1); else
            interpAvg.trial = interpAvg.trial(:,:,1); end

        %     figure, subplot(121), plot(1:sz(2), squeeze(data.trial(idx(1),:,1))); % visualise trial
        %     subplot(122), plot(1:sz(2), squeeze(interpAvg.avg ) ); % visualise interpolated trial
        %     ylim([-20 30])
        % sanity check that the labels match
        lv_match_channels(data, interpAvg.label)

        if isfield(interpAvg,'avg'), data_interp.trial(idx,:,1) = interpAvg.avg; inter_co=inter_co+length(idx); else
            data_interp.trial(idx,:,1) = interpAvg.trial; inter_co=inter_co+length(idx); end
    end
    toc

    % reshape back
    data_interp.trial = data_interp.trial(:,:,1);
    data_interp.trial = reshape(data_interp.trial',sz(2),sz(1),sz(3));
    data_interp.trial = permute(data_interp.trial,[2 1 3]);
    data_interp.time = temp_time;
elseif choice==0 % now the case where we get the continuous signal of every channel and see
    % if it should be removed and also the comparison between its variance
    % and the variance of other channels
    THRESHOLD = 5;
    data_interp = data; lv_layout = data.lv_layout;
    % variance channels vs. each other
    dat = permute(data.trial,[2 3 1]); dat = reshape(dat,size(data.trial,2),[]); %2d
    all_measure = var(dat,[],2);
    iqr_val = iqr(all_measure); q1=prctile(all_measure,25); q3=prctile(all_measure,75);
    suggested_bad_channel = data.label(all_measure<q1-(THRESHOLD*iqr_val) | all_measure>q3+(THRESHOLD*iqr_val));
    data_interp.suggested_bad_channel=suggested_bad_channel;
    % continuous measure of outliers within the same channel
    for i=1:size(dat,1)
        all_measure = dat(i,:);
        iqr_val = iqr(all_measure); q1=prctile(all_measure,25); q3=prctile(all_measure,75);
        bad_flag(i,:) = double(all_measure<q1-(THRESHOLD*iqr_val) | all_measure>q3+(THRESHOLD*iqr_val));
    end
    bad_flag = reshape(bad_flag,size(data.trial,2),size(data.trial,3),[]); bad_flag = permute(bad_flag,[3 1 2]);
    bad_flag = sum(bad_flag,3); goodFlags_matrix = ~bad_flag;
    % visualise bad segments
    idx = find(sum(bad_flag,2)==0);
    figure, subplot(121), plot(0:1/200:(size(dat,2)-1)/200, dat); title('all trials');
    subplot(122), plot(0:1/200:(size(dat,2)-1)/200, dat); title('all trials');
    temp=data.trial; temp(idx,:,:)=nan;
    temp = permute(temp,[2 3 1]); temp = reshape(temp,size(data.trial,2),[]); %2d
    hold on, plot(0:1/200:(size(temp,2)-1)/200, temp, 'r');

    disp(suggested_bad_channel);

    temp=dat;
    temp(~ismember(data.label,data_interp.suggested_bad_channel),:)=nan;
    hold on, plot(0:1/200:(size(dat,2)-1)/200, temp, 'k'); %ylim([-400 400]);
    % fix the bad channels?
    suggested_bad_channel=[];% to force to not correct the bad channel and we will remove it later..
    if ~isempty(suggested_bad_channel)
        if strcmp(data.fix_option,'auto')==1 opt = 1; warning('accepting channel fix..'); else 
        opt = listdlg('ListString',{'yes','no','visualinspection'},'ListSize',[450,450]); end
        if opt==1 % then this black channels are bad
            goodFlags_matrix(:,ismember(data.label,data_interp.suggested_bad_channel),:)=0;
        elseif opt==3
            cfg=[]; cfg.continuous ='yes';  cfg.ylim = [-400 400];
            cfg.blocksize=size(data.trial,1)*data.time(end); % check if the idenitified bad channel is the same as the visual
            ft_databrowser(cfg,data);
            ch = lv_tune_params('Channel','');
            goodFlags_matrix(:,ismember(data.label,ch),:)=0;
        end
    end

    % correct bad channels
    % Find neighbours
    cfg              = [];
    cfg.method       = 'triangulation'; % the way it's going to connect channels indicating neighborhood
    cfg.senstype     = 'EEG';
    cfg.layout       = lv_layout; cfg.feedback     = 'no';
    neighbours	= ft_prepare_neighbours(cfg);
    tic
    [unique_trls,~,ic] = unique(goodFlags_matrix,'rows');
    inter_co = 0; % interpolation count to make sure it was done correctly
    data_interp.suggested_bad_trial = [];
    for i=1:size(unique_trls,1)
        idx = find(ic == i);

        cfg                  = [];
        cfg.badchannel       = data.label( ~goodFlags_matrix(idx(1),:) );
        if isempty(cfg.badchannel), data_interp.trial(idx,:,:) = data.trial(idx,:,:); continue; end % good trial
        if (length(cfg.badchannel)>0.25*length(data.label)), data_interp.suggested_bad_trial=[data_interp.suggested_bad_trial idx(:)']; end
        cfg.neighbourdist    = 4;
        cfg.neighbours     	 = neighbours;
        cfg.layout           = lv_layout;
        cfg.method           = 'spline';
        cfg.feedback         = 'no';
        cfg.trackcallinfo	 = 'no';
        cfg.trials           = idx;

        interpAvg = ft_channelrepair(cfg ,data);

        %     figure, subplot(121), plot(1:sz(2), squeeze(data.trial(idx(1),:,1))); % visualise trial
        %     subplot(122), plot(1:sz(2), squeeze(interpAvg.avg ) ); % visualise interpolated trial
        %     ylim([-20 30])
        % sanity check that the labels match
        lv_match_channels(data, interpAvg.label)

        if isfield(interpAvg,'avg'), data_interp.trial(idx,:,:) = interpAvg.avg; inter_co=inter_co+length(idx); else
            data_interp.trial(idx,:,:) = interpAvg.trial; inter_co=inter_co+length(idx); end
    end
    toc

    % visualising data
    dat = permute(data_interp.trial,[2 3 1]); dat = reshape(dat,size(data_interp.trial,2),[]); %2d
    figure, plot(0:1/200:(size(dat,2)-1)/200, dat); title('all trials');
    ylim([-400 400]);

    if strcmp(data.fix_option,'auto')==1 close all; end
end



end



function  lv_match_channels(segmented_data, label)
% check that the labels are exactly the same
newIdx=[];
for i=1:length(label)
    newIdx = [newIdx ; find( ismember(segmented_data.label , label{i}) )];
end

idx = 1:length(label);
if any( idx' - newIdx  )
    error('Channel mismatch');
end
end


