function [ varargout ] = lv_post_classification(cfg)
% takes cfg and the cfg should even contain the data or
% anything that might be needed according to the cfg.method used ...

tested =1; % tested 2d only
dbstop in lv_post_classification at 10 if tested==0


method = cfg.method;

fprintf(['\n Performing post classification: ' method ' \n']);

%%
switch method

    case 'blocks_performance' % 1D or 2D blocks of trials performance through out stimulation time
        %%
        % signi_mask should contain one cluster at a time .. and takes trn
        % and tst data sets, returns a performance curve across stimulation time
        [dimensions,blocks_perc,shift_perc,slide_trn] = lv_tune_params('do you want 1D classification?','0','percentage of trials per block','20',...
            'shifting between blocks, percentage','1','are we sliding on train data?','1', 'just use defaults');


        data_trn = cfg.data_trn; data_tst = cfg.data_tst;


        warning('be careful assumed classifier is: mvpa,timextime classification');
        if slide_trn==1, sz = size(data_trn.trial); else,  sz = size(data_tst.trial); end % trl_ch_time
        signi_mask = cfg.signi_mask; % 2d of binary mask
        signi_mask(signi_mask~=0)=1; signi_mask(signi_mask==0)=nan;
        [rectangle_row,rectangle_col] = find(signi_mask==1); % find the rectangle containing the reactivation and aplly classification only inside it
        rectangle_row = min(rectangle_row):max(rectangle_row);
        rectangle_col = min(rectangle_col):max(rectangle_col);

        data_trn.trial = data_trn.trial(:,:,rectangle_row);
        data_tst.trial = data_tst.trial(:,:,rectangle_col);

        blocks = ceil(sz(1)*(blocks_perc/100)); shift = ceil(sz(1)*(shift_perc/100));
        trials_to_pick = (1:shift:(sz(1)-blocks))'; trials_to_pick = [trials_to_pick trials_to_pick+blocks]; %2d of every row is the trials to pick

        %         temp = round(length(data_trn.trial)/2);
        %         trials_to_pick = [1 temp; temp+1 length(data_trn.trial)];

        trials_to_pick = [1 length(data_trn.trial)];

        cfg_classifier = [];
        cfg_classifier.classifier      = 'lda';
        cfg_classifier.metric          = 'auc'; % 'dval' distance from boundary 23taked
        cfg_classifier.preprocess      = 'undersample';
        cfg_classifier.hyperparameter           = [];
        cfg_classifier.hyperparameter.lambda    = 'auto';

        for i=1:size(trials_to_pick,1)
            acc{1,i} = nan(size(signi_mask));
            keep = trials_to_pick(i,1):trials_to_pick(i,2);
            if slide_trn==1 % taking subset of trials according to the set to be segmented into blocks
                training = data_trn.trial(keep,:,:); training_lbl = data_trn.trialinfo(keep);
                testing = data_tst.trial; testing_lbl = data_tst.trialinfo;
            else
                training = data_trn.trial; training_lbl = cfg.trialinfo;
                testing = data_tst.trial(keep,:,:); testing_lbl = data_tst.trialinfo(keep);
            end

            rng(1)
            if dimensions==0 % 2D
                warning('be careful assumed classifier is: mvpa,timextime classification (lda,acc,undersample), with no z-scoring you should z-score beforehand');
                [acc{1,i}(rectangle_row,rectangle_col),~ ] = mv_classify_timextime(cfg_classifier,...
                    training , training_lbl , testing  , testing_lbl); % performance is put in the rectangle location only
            else             % 1D
                warning('be careful assumed classifier is: mvpa,classify_across_time classification (lda,acc,undersample), with no z-scoring you should z-score beforehand');
                [acc{1,i}(rectangle_row,rectangle_col),~ ] = mv_classify_across_time(cfg_classifier,...
                    training , training_lbl , testing  , testing_lbl);
            end
            acc{1,i} = acc{1,i} .* signi_mask;
            blocks_curve(1,i) = nanmean(acc{1,i}(:));
        end

        % normalizing the time to 0->1 .. because the trials lengths are different for diff. sbjs
        %         pts = linspace(0,1,length(blocks_curve)); resolution=0.01; warning(['resolution is ' num2str(resolution) ' use it for plotting 0:resolution:1']);
        %         normalized_block_curve = interp1(pts, blocks_curve  , 0:resolution:1);

        varargout{1} = blocks_curve;


    case 'reactivations_types'
        %% if wake is trn and sleep is test we can directly extract the index of the correct trials of sleep and see if they reoccurring acc. to the time of the peaks
        % but if sleep is the training then we a method similar to RSA so
        % we take every trial x and see its class and see if it correlates
        % signi. with the majority of wake trials that's how we say that
        % it's pure .. we do that for x at the time of the 1st peak and
        % then repeat at the second peak and see how many trials with just
        % pk1 .. just pk2 .. reoccurring..
        if cfg.slp_trn==1
            data = cfg.sleep; sleep_lbls = cfg.sleep.trialinfo;
            im = cfg.im;  im_lbls = cfg.im.trialinfo;
            signi_mask = cfg.signi_mask; % 2d of binary mask

            min_size_cluster = 1:size(signi_mask,1); % assume min. size cluster and then update it with the new min.
            L = bwlabel(signi_mask,8); % find the connected points in any direction '8-connectivity' and labels the connected clusters

            trials_types_flags = zeros(length(sleep_lbls),1);
            for i=1:(max(L(:)))
                [row_col{i,1},row_col{i,2}] = find(L==i); % row_col: clusters_rows,and,columns
                row_col{i,1} = unique(row_col{i,1});  row_col{i,2} = unique(row_col{i,2});

                wake_trl = squeeze(im.trial(:,:,row_col{i,2}));
                % reducing duration of wake
                pts = linspace(1,size(wake_trl,3),length(row_col{i,1}) );
                compressed_wake = [];
                for trl=1:size(wake_trl,1)
                    for ch=1:size(wake_trl,2)
                        compressed_wake(trl,ch,:) = interp1(1:size(wake_trl,3), (squeeze(wake_trl(trl,ch,:)))' ,pts);
                    end
                end

                sleep_dat = squeeze(data.trial(:,:,row_col{i,1})); sleep_dat = permute(sleep_dat, [3 2 1]); % time_ch_trl
                sleep_dat = reshape(sleep_dat, size(sleep_dat,1)*size(sleep_dat,2),[] ); % TimeCh_trl

                im_dat = permute(compressed_wake, [3 2 1]); % time_ch_trl
                im_dat = reshape(im_dat, size(im_dat,1)*size(im_dat,2),[] ); % TimeCh_trl


                [rho, pval] = corr(sleep_dat, im_dat, 'type','spearman'); % ex: 1st row the frst trl from sleep corr with all trls from wake

                % just to visualize
                class1_rho = rho(sleep_lbls==1,im_lbls==1);  class2_rho = rho(sleep_lbls==2,im_lbls==2);


                % pk1 .. trials of sleep correlating signi. with wake for
                % its correct class .. peak1 is any peak according to the loop so I just named it pk1 for simplicity
                % left hand
                true_class1_trls_idx = find(sleep_lbls==1);
                % majority correlating signi. and positively with the correct class .. was stringent so now it's: more trls with correct than incorrect
                [pk1_class1_sleep,~] = find ( sum( rho(sleep_lbls==1,im_lbls==1)>0 & pval(sleep_lbls==1,im_lbls==1)<0.05 , 2) ...
                    > sum( rho(sleep_lbls==1,im_lbls==2)>0 & pval(sleep_lbls==1,im_lbls==2)<0.05 , 2)  );% majority: (size(rho,2)/2) );
                pk1_class1_sleep = true_class1_trls_idx(pk1_class1_sleep);

                % right hand
                true_class2_trls_idx = find(sleep_lbls==2);
                % majority correlating signi. and positively with the correct class .. was stringent so now it's: more trls with correct than incorrect
                [pk1_class2_sleep,~] = find ( sum( rho(sleep_lbls==2,im_lbls==2)>0 & pval(sleep_lbls==2,im_lbls==2)<0.05 , 2) ...
                    > sum( rho(sleep_lbls==2,im_lbls==1)>0 & pval(sleep_lbls==2,im_lbls==1)<0.05 , 2)  );% majority: (size(rho,2)/2) );
                pk1_class2_sleep = true_class2_trls_idx(pk1_class2_sleep);

                flag = i; % pk1 will be 1 and pk2 will be 2
                trials_types_flags(pk1_class1_sleep) = trials_types_flags(pk1_class1_sleep) + flag; % such that when it includes both it will be 3
                trials_types_flags(pk1_class2_sleep) = trials_types_flags(pk1_class2_sleep) + flag;

            end

            varargout{1} = trials_types_flags; % matrix (trls_clusters) of flags with pk1:1 .. pk2:2 .. reoccurring:3 ,, non:0
        end
    case 'are the reactivations happenning in the same trial when training is with sleep?'
        %% since training is with sleep we want to see inside the trial itself and correlate the significant clusters together
        % % reduce the length of the big cluster(s) to match the smallest one and correlate with the smallest one .. to
        % balance the number of samples in the correlation

        % this is a stringent test because it needs the trials to correlate
        % at the two peaks and there is no guarantee for that but if it
        % worked then they are truly reoccurring but if it didn't then you
        % cannot say for sure that they aren't reoccurring .. we didn't use
        % it eventually because of being very stringent and used
        % 'reactivations_types' to check the number of significant and
        % positive correlation if being correct more for the correct label
        % rather than being correct for the majority of trials (this one
        % is sleep sleep because it needs to correlate the peaks with each
        % other)

        data = cfg.data;
        signi_mask = cfg.signi_mask; % 2d of binary mask

        min_size_cluster = 1:size(signi_mask,1); % assume min. size cluster and then update it with the new min.
        L = bwlabel(signi_mask,8); % find the connected points in any direction '8-connectivity' and labels the connected clusters
        for i=1:(max(L(:)))
            [row_col{i,1},row_col{i,2}] = find(L==i); % row_col: clusters_rows,and,columns
            row_col{i,1} = unique(row_col{i,1});  row_col{i,2} = unique(row_col{i,2});
            if length(row_col{i,1})<length(min_size_cluster), min_size_cluster = row_col{i,1}; clusternum=i; end % min_size_cluster is row length
        end


        for i=1:(max(L(:)))
            if i==clusternum, continue; end

            trials_with_reoccur(i,1)=0;
            for trl=1:size(data.trial,1) % correlate the trials at the time of clusters that we got
                for ch=1:size(data.trial,2)
                    x = squeeze(data.trial(trl,ch,row_col{i,1}));
                    pts = linspace(1,length(x),length(min_size_cluster));
                    compressed_x = interp1(1:length(x), x ,pts);
                    [correlations(i,trl,ch,1),correlations(i,trl,ch,2)] = corr(squeeze(data.trial(trl,ch,min_size_cluster)),compressed_x','type','spearman');
                    % cluster_trial_channel_r cluster_trial_channel_p
                end
                if sum(correlations(i,trl,:,2)<0.05 & correlations(i,trl,:,1)>0) > floor(size(data.trial,2)/2)
                    % count the channels with significant similarity between the clusters if more than half then
                    % consider the trial as having the same effect
                    trials_with_reoccur(i,1) = trials_with_reoccur(i,1)+1;
                end
            end

        end

        varargout{1} = (trials_with_reoccur / size(data.trial,1))*100; % percentage of trials with reoccurence...

    case 'recurrence plot'
        if strcmp(cfg.approach,'correlation')~=1
            % this will return a recurrence plot that should reveal
            % recurrence and it's based on the values when they are almost
            % the same this will reflect higher recurrence and it's done in
            % a timextime manner

            % but before using it ask is what is recurring reactivation?
            % because there could be recurrence because the pattern is
            % recurring but doesn't indicate reactivation it could erp
            % pattern that's recurring..
            threshold = 5; % value for recurrence consideration (if <=threshold apart then they are considered the same)
            data = cfg.data;
            for i=1:size(data.trial,1)
                for j=1:size(data.trial,2)
                    vec = squeeze(data(i,j,:));
                    for tt=1:size(data.trial,3)
                        recurr(i,j,tt,:) = (vec - vec(tt)) <= threshold; % trl_ch_time_time
                        recurr(i,j,tt,:) = norm(vec - vec(tt)) <= threshold;
                    end
                end
            end

            varargout{1} = recurr; % trl_ch_time_time .. recurrence plot
        else
            data = cfg.data;
            % based on correlation
            agg_similarity = zeros(size(data.trial,3), size(data.trial,3)); % agregating similarities from all trials
            for i=1:size(data.trial,1)
                agg_similarity = agg_similarity + corr(squeeze(data.trial(i,:,:)), squeeze(data.trial(i,:,:)));
            end
            varargout{1} = agg_similarity./size(data.trial,1);
        end

    case 'compression_correct_pts'
        % when we have a matrix with trialxtime containing flags with the
        % correct pts after classification and we want to compare the
        % length between the correct pts in wake and sleep and fill the
        % small gaps <40ms between trust worthy reactivations >40ms
        lens_sleep=[]; recurrence_distrib=[];
        sleep_res = cfg.data;
        for i=1:size(sleep_res,1)
            v = sleep_res(i,:);
%             % fixing gaps .. smoothing fixes this automatically because it considers the majority
%             CC = bwconncomp(v);
%             temp = cell2mat(cellfun(@(x) [x(1);x(end)], CC.PixelIdxList,'Un',0));
%             temp2 = temp(:); temp2=diff(temp2); id=1:length(temp2); id = mod(id,2);
%             good_lens = temp2(id==1)+1; gaps = temp2(id==0)-1; % not inclusive of edges
% 
%             tofill = gaps <= 5; % 8*5=40ms
%             good_lens = conv(good_lens',[1 1],'valid') >= 10; % summing each two consecutive elements
%             tofill = find(tofill(:) & good_lens(:)); % lengths are >40ms and the gap is <40ms
%             for j=1:length(tofill), v( temp(2,tofill(j))+1 : temp(1,tofill(j)+1)-1 ) = 1; end

            CC = bwconncomp(v); temp = cell2mat(cellfun(@(x) length(x), CC.PixelIdxList,'Un',0));
            lens_sleep = [lens_sleep temp];
            if cfg.recurrence==1
                recurrence_distrib = [recurrence_distrib sum(ismember(temp,cfg.ratios))];
            end
        end
        varargout{1} = lens_sleep;
        if cfg.recurrence==1, varargout{2}=recurrence_distrib; else varargout{2}=[]; end

    case 'compression_wilson'
        % Wislon's replay compression by varying the size of the long
        % dataset and performing rsa/classification and then looking at the
        % peak of classification to determine the compression ratio
        % Wilson's used a scaling factor (SF) which is the denominator of the division of lengths (length short/SF)
        % so because encoding/run in his data was long he divided the
        % length of rem by the SF so a SF=2 means that replay is dilated in
        % sleep because length of sleep/2 is very small so he gets this
        % small length to be the length of the window of wake and then
        % resizes this to match the length of sleep episode, so this means it was compressed in the run so dilated in sleep
        data_long = cfg.data_long; % here long data is sleep for mri data which is the opposite of wilson's
        data_short = cfg.data_short;
        centered = cfg.centered;
        SF_samples = round(logspace(log10(10),log10(length(data_short.time)),10)); SF_samples(end)=[]; % compression .. log scaling will give more values
        % near each other when we do division because we want to look at high compressions and not have huge jumps
        SF_samples = [SF_samples round(logspace(log10(length(data_short.time)),log10(length(data_long.time)),10))]; % dilation
        tic
        for i=1:length(SF_samples)
            cfg=[]; cfg.data = data_long;
            cfg.shift = 1; % shifts in samples between windows (minimum is one sample shift because then we will get the full resolution)
            cfg.windlen = SF_samples(i); % because we take the samples around the current sample, 
            % so the actual window is wider so we consider this later at varargout{1}.compressionRatio
            result = lv_slider(cfg); 
            
            progressbar = ParforProgressbar(size(result.higher_data,4),'title', ['compression ' num2str(i) ' of ' num2str(length(SF_samples))]);
            sz = size(data_short.trial); data_short_hold=data_short;
            parfor j=1:size(result.higher_data,4) % looping on sliding windows
                temp = squeeze(result.higher_data(:,:,:,j));
                window = imresize3(temp,[size(temp,1) size(temp,2) sz(3)]); % resize each chxtime
                 
                % compressing ch_time with PCA
                cfg=[]; cfg.data.trial = reshape(window,size(window,1),[]); temp2=[];
                temp2(1,:,:)=(cfg.data.trial)'; 
                cfg.data.trial=[]; cfg.data.trial=temp2; % because data is 2d we simulate trl_ch_time
                cfg.method = 'pca'; cfg.step = 'calculate'; cfg.centered=centered;
                comp = lv_component_analysis(cfg); temp2=[];
                % transform
                cfg.eigVects = comp.eigVects; cfg.chosen = zeros(length(comp.eigVals),1);
                id = find(cumsum(comp.eigVals)> 0.95);
                cfg.chosen(1:id(1))=1;
                cfg.step = 'transform';
                window = lv_component_analysis(cfg);
                window = squeeze(window)'; if sum(cfg.chosen)==1, window=window(:); end
                cfg.data.trial = reshape(data_short_hold.trial,size(data_short_hold.trial,1),[]);
                temp2(1,:,:)=(cfg.data.trial)'; 
                cfg.data.trial=[]; cfg.data.trial=temp2; temp2=[];
                temp2 = lv_component_analysis(cfg);
                temp2 = squeeze(temp2)'; if sum(cfg.chosen)==1, temp2=temp2(:); end
                
                % Classification
                cfg=[]; cfg.method = 'timextime';
                cfg.classifier_type = {'lda'}; cfg.perf_measure='acc';
                cfg.trn.trialinfo=data_long.trialinfo(:,1);
                cfg.tst.trial = temp2; cfg.tst.trialinfo=data_short.trialinfo(:,1);
                cfg.folds=nan; cfg.do_parallel=0;

                cfg.trn.trial = window; 
                cfg.trn.trial = zscore(cfg.trn.trial,[],1);  cfg.tst.trial = zscore(cfg.tst.trial,[],1); 
                classification_pk(i,j) = lv_classify(cfg); % SFs x sliding_windows 
                temp2=[];
                progressbar.increment();    
            end 
            delete(progressbar);
        end
        toc
        varargout{1}.score = max(classification_pk,[],2); % the peak of all sliding windows is the score for that SF
        varargout{1}.compressionRatio = (SF_samples+1)./length(data_short.time); % e.g., if the peak at 2 samples means that two samples from the long data
        % when resized to the short data length is peaking .. so
        % reactivation is compressed in the long data and the compressionRatio is 2/length_short_data
        % or swap it length_short_data/2 to get how much it's faster

    case 'compression_wilson_jittered'
        % similar to the previous case but in this one we put all sliding
        % windows as trials because we don't know when the effect happens
        % so we can't time lock trials together in the sliding window ..
        % also considers long trials with nans and then removes missing values
        data_long = cfg.data_long; % here long data is sleep for mri data which is the opposite of wilson's
        data_short = cfg.data_short; 
        SF_samples = round(logspace(log10(10),log10(length(data_short.time)),10)); % compression .. log scaling will give more values
        % near each other when we do division because we want to look at high compressions and not have huge jumps
        SF_samples = [SF_samples round(length(data_short.time) + cumsum(flip( diff(logspace(log10(length(data_short.time)),log10(length(data_long.time)),10)) )))];
 
        tic
        for i=1:length(SF_samples)
            cfg=[]; cfg.data = data_long;
            cfg.shift = 2; % shifts in samples between windows (minimum is one sample shift because then we will get the full resolution)
            cfg.windlen = SF_samples(i); % because we take the samples around the current sample, 
            % so the actual window is wider so we consider this later at varargout{1}.compressionRatio
            result = lv_slider(cfg); % trl_ch_time_windows 
            sz = size(data_short.trial); 
            temp = permute(result.higher_data,[2 3 1 4]); % ch_time_trl_windows
            temp = reshape(temp,size(temp,1),size(temp,2),[]);
            window = permute(temp, [3 1 2]);
            window = imresize3(window,[size(window,1) size(window,2) sz(3)]);
            clabels = repmat(data_long.trialinfo(:,1),size(result.higher_data,4),1); 
            
            % PCA 
            comp=[]; window = reshape(window,size(window,1),[]);
            id = find(isnan(sum(window,2))); % removing nans
            window(id,:)=[]; clabels(id)=[];
            window = window - repmat(mean(window,1),size(window,1),1); % centering
            data_short_mat = reshape(data_short.trial,size(data_short.trial,1),[]);
            data_short_mat = data_short_mat - repmat(mean(data_short_mat,1),size(data_short_mat,1),1); % centering
            [comp.eigVects,~,comp.eigVals] = pca(window); % gets eig vectors and values and they are sorted
            comp.eigVals = comp.eigVals./sum(comp.eigVals); 
            id = find(cumsum(comp.eigVals)> 0.95);
            window = window * comp.eigVects(:,1:id(1)); % transform window
            data_short_mat = data_short_mat * comp.eigVects(:,1:id(1)); 

            % Classification
            cfg=[]; cfg.method = 'timextime';
            cfg.classifier_type = {'lda'}; cfg.perf_measure='acc';
            cfg.trn.trialinfo= clabels;
            cfg.tst.trial = data_short_mat; cfg.tst.trialinfo=data_short.trialinfo(:,1);
            cfg.folds=nan; cfg.do_parallel=1;

            cfg.trn.trial = window;
            cfg.trn.trial = zscore(cfg.trn.trial,[],1);  cfg.tst.trial = zscore(cfg.tst.trial,[],1);
            classification_pk(i,1) = lv_classify(cfg); % SFs x sliding_windows
            lv_progress(i,length(SF_samples),'Compression/dilation: ');
        end
        toc
        varargout{1}.score = classification_pk;
        varargout{1}.compressionRatio = (SF_samples+1)./length(data_short.time); % e.g., if the peak at 2 samples means that two samples from the long data
        % when resized to the short data length is peaking .. so
        % reactivation is compressed in the long data and the compressionRatio is 2/length_short_data
        % or swap it length_short_data/2 to get how much it's faster

end








