function data = lv_align(cfg)
% aligns signals inside .trial to have the maximum correlation values
% uses reference from the same data and the reference is chosen in the main
% function code ... EXAMPLE at the end of this code.
data = cfg.data; % data.time all time extending beyond trial limit
roi = cfg.roi; % the trial time i.e.,: [0 1.1]
data.trial = (data.trial);
maxShift = (cfg.shift*200)/1000;% 50ms both directions so 100ms proximity
% if peaks are 100ms apart they shall still be aligned
roi_id = [nearest(data.time,roi(1)) nearest(data.time,roi(2))];

% shifting channels together is the correct way because otherwise we will
% make new signals that don't exist if we shift each channel individually
% and also every channel won't have the same temporal alignment of the
% other channels because it is optimised within itself so when the
% classifier works on timepts to classify it won't get consistent good features..
% with matlab's normalised correlation which is 2d and can do 1d as in the example below
% example ... this example was edited to be 1D but compare it to the
% example inside normxcorr2 and the result should be similar
%         peppers = im2gray(imread('peppers.png'));
%         onion = peppers(:,140:280);
%         montage({peppers,onion})
%
%         c = normxcorr2(onion,peppers);
%         surf(c)
%         shading flat
%
%         % only looking at the logical parts without edge artifacts of the 2D convolution, so we remove the padded convolution results
%         c(:,[1:size(onion,2)-1 end-(size(onion,2)-1)+1:end])=[];
%         c([1:size(onion,1)-1 end-(size(onion,1)-1)+1:end],:)=[];
%
%         [~,xpeak] = max(c);
%
%         imagesc(peppers)
%         drawrectangle(gca,'Position',[xpeak,1,size(onion,2),size(onion,1)], ...
%             'FaceAlpha',0);

% trls = data.trial(:,:,roi_id(1):roi_id(2)); sz = size(trls);
longer_trls = data.trial(:,:,roi_id(1)-maxShift:roi_id(2)+maxShift);
lag = [];
sz = size(longer_trls);
% progressbar = ParforProgressbar(size(data.trial,1),'title', 'trials'' correlations progress:');
% parfor i=1:size(data.trial,1)
for j=1:size(data.trial,1)
%     c = normxcorr2( squeeze(trls(i,:,:)), squeeze(longer_trls(j,:,:)) ); % j is the longer and the one that enables shifting
    c = normxcorr2( cfg.ref_trl, reshape(squeeze(longer_trls(j,:,:)),sz(2),sz(3)) );
    % only looking at the logical parts without edge artifacts of the 2D convolution, so we remove the padded convolution results
    c(:,[1:size(cfg.ref_trl,2)-1 end-(size(cfg.ref_trl,2)-1)+1:end])=[];
    c([1:size(cfg.ref_trl,1)-1 end-(size(cfg.ref_trl,1)-1)+1:end],:)=[];
    if isfield(cfg,'classful'), [~,xpeak] = max(c); else, [~,xpeak] = max(abs(c)); end % abs for classless correlation,  classful for same shape align like for group level align
    
    %lag(j,i) = -(xpeak - maxShift -1); % will reveal the shift amount by -maxShift because if the trl won't be shifted it will appear here
    % as shifted by maxShift because we appended the longer_trls... it
    % is negative to align them correctly
    lag = [lag   -(xpeak - maxShift -1)]; % see the commented part above we just do it in a row for parfor to work
    %        lag(i,j) = -lag(j,i);
end
%     progressbar.increment();
% end
% delete(progressbar);

% lags=zeros(size(data.trial,1) ,size(data.trial,1)); count=0; % refilling with the results from parfor
% for i=1:size(data.trial,1)
%     for j=i:size(data.trial,1)
%         count=count+1; lags(j,i) = lag(count);
%         lags(i,j) = -lags(j,i);
%     end
% end
% 
% 
% lag = round(mean(lags,2));
lag=lag(:);
if any(lag>maxShift), error('shifting by more than the maxShift !'); end
% shifting trials according to the result
for i=1:size(data.trial,1)
    lv_progress(i,size(data.trial,1),'aligning progress: ');
    for ch=1:size(data.trial,2)
        if size(lag,2)==1, lag = repmat(lag,1,size(data.trial,2));  end % repeating because that's topos so all channels will be shifted together
        data.trial(i,ch,:) = circshift(data.trial(i,ch,:),lag(i,ch));
    end
end
% data.trial = data.trial(:,:,roi_id(1):roi_id(2));
% data.time = data.time(roi_id(1):roi_id(2));

data.trial(:,:,1:maxShift)=[]; data.trial(:,:,end-(maxShift-1):end)=[]; % removing the parts with circular shift
data.time(1:maxShift)=[]; data.time(end-(maxShift-1):end)=[];
%% other ways


% lags = -maxShift:maxShift;
% progressbar = ParforProgressbar(size(data.trial,1),'title', 'trials'' correlations progress:');
% parfor trl=1:size(data.trial,1)
%     similarity=[]; %lv_progress(trl,size(data.trial,1),'trials progress: ');
%     all_trls = data.trial(:,:,roi_id(1):roi_id(2));
%     for i=1:length(lags)
%         if ~isfield(cfg,'method') % lag is inverted because when we move with ids we get the right most so to put the right most in the left part that will be a negative circular shift..
%             similarity(i,:,:) = sum(bsxfun(@times, data.trial(trl,:,roi_id(1)-lags(i):roi_id(2)-lags(i)), all_trls) ,3);
%         else % lag_alltrls_ch
%             if strcmp(cfg.method,'topos'), similarity(i,:,:) = sum(sum(bsxfun(@times, data.trial(trl,:,roi_id(1)-lags(i):roi_id(2)-lags(i)), all_trls) ,3),2);
%             end, end % we compress on channels as well so a whole topo not individual channels
%     end
%     [~,best_lag] = max(abs(similarity),[],1); % alltrials_ch
%     if isfield(cfg,'method'), if strcmp(cfg.method,'topos'), best_lag=best_lag'; end, end
%     apply_shift(trl,:) = lags(round(mean(squeeze(best_lag),1))); %  1_ch .. mean best lag for every channel ch
%     progressbar.increment();
% end
% delete(progressbar);
% % shifting trials according to the result
% for i=1:size(data.trial,1)
%     lv_progress(i,size(data.trial,1),'aligning progress: ');
%     for ch=1:size(data.trial,2)
%         if size(apply_shift,2)==1, apply_shift = repmat(apply_shift,1,size(data.trial,2));  end % repeating because that's topos so all channels will be shifted together
%         data.trial(i,ch,:) = circshift(data.trial(i,ch,:),apply_shift(i,ch));
%     end
% end
% data.trial = data.trial(:,:,roi_id(1):roi_id(2));
% data.time = data.time(roi_id(1):roi_id(2));
% end



% with different ways to pick the fast one
%  for trl=1:size(data.trial,1)
%     similarity=[]; lv_progress(trl,size(data.trial,1),'trials progress: ');
%     for i=1:length(lags)
% %         trial_shifted = repmat(  data.trial(trl,:,roi_id(1)-lags(i):roi_id(2)-lags(i)), size(data.trial,1),1,1 ); % shifted ch_time slice.. -lag because the shift that we do is in the other direction
% %         similarity(i,:,:) = sum((trial_shifted .* data.trial(:,:,roi_id(1):roi_id(2))),3); % lag_alltrials_ch
%         similarity(i,:,:) = sum(bsxfun(@times, data.trial(trl,:,roi_id(1)-lags(i):roi_id(2)-lags(i)), data.trial(:,:,roi_id(1):roi_id(2))) ,3);
%     end
% %     for i=1:length(lags)
% %         trial_shifted = squeeze(data.trial(trl,:,roi_id(1)-lags(i):roi_id(2)-lags(i))); % shifted ch_time slice.. -lag because the shift that we do is in the other direction
% %         for ch=1:size(trial_shifted,1)
% %             similarity(i,:,ch) = trial_shifted(ch,:) * squeeze(data.trial(:,ch,roi_id(1):roi_id(2)))'; % lag_trl_alltrials_ch
% %         end
% %     end
%     [~,best_lag] = max(abs(similarity),[],1); % alltrials_ch
%     apply_shift(trl,:) = lags(round(mean(squeeze(best_lag),1))); %  1_ch .. mean best lag for every channel ch
% end



%% EXAMPLE
%% trial level
%     % lv_align choosing a reference signal (classless it should be a global reference
%     toi_len=diff(train_ival)/2;
%     cfg=[]; cfg.latency = [ 0-toi_len 1.5+toi_len]; data = ft_selectdata(cfg,sleep); coeff=[];
%     progressbar = ParforProgressbar(size(data.trial,1),'title', 'trials'' correlations progress:');
%     parfor i=1:size(data.trial,1)
%         for j=i:size(data.trial,1)
%             c = normxcorr2( squeeze(data.trial(i,:,:)), squeeze(data.trial(j,:,:)) );
%             
%             % only looking at the logical parts without edge artifacts of the 2D convolution, so we remove the padded convolution results
%             c(:,[1:size(data.trial,3)-1 end-(size(data.trial,3)-1)+1:end])=[];
%             c([1:size(data.trial,2)-1 end-(size(data.trial,2)-1)+1:end],:)=[];
%             coeff = [coeff   abs(c)]; % see the commented part above we just do it in a row for parfor to work
%             %        lag(i,j) = -lag(j,i);
%         end
%         progressbar.increment();
%     end
%     delete(progressbar);
%     count=0; % refilling with the results from parfor
%     coeff_all=[];
%     for i=1:size(data.trial,1)
%         for j=i:size(data.trial,1)
%             count=count+1; coeff_all(j,i) = coeff(count);
%             coeff_all(i,j) = coeff_all(j,i);
%         end
%     end
%     
%     coeff = mean(coeff_all,2);
%     
%     [val,ref_trl(nn,1)] = max(coeff)
%     % aligning with lv with xcorrelation
%     cfg=[];
%     cfg.data = sleep; % data.time all time extending beyond trial limit
%     cfg.roi =[0-toi_len 1.5+toi_len]; cfg.method='topos';
%     cfg.shift = 200;
%     cfg.ref_trl = squeeze(data.trial(ref_trl(nn,1),:,:));
%     sleep = lv_align(cfg);
 %% lv_align sbj level classful align (same shape align)
% % choosing a reference signal 
% auc = (squeeze(tempAccM))';
% toi_len=diff(train_ival)/2; 
% ids = nearest(xax, 0-toi_len ):nearest(xax, 1.5+toi_len );
% data=[]; data.trial(:,1,:) = auc(:,ids); data.time = xax(ids); 
% 
% coeff=[]; 
% for i=1:size(data.trial,1)
%     for j=i:size(data.trial,1)
%         c = normxcorr2( squeeze(data.trial(i,:,:)), squeeze(data.trial(j,:,:)) );   c=c';     
%         % only looking at the logical parts without edge artifacts of the 2D convolution, so we remove the padded convolution results
%         c(:,[1:size(data.trial,3)-1 end-(size(data.trial,3)-1)+1:end])=[]; 
%         
%         coeff = [coeff   abs(c)]; % see the commented part above we just do it in a row for parfor to work
%         %        lag(i,j) = -lag(j,i);
%     end 
% end 
% count=0; % refilling with the results from parfor
% coeff_all=[];
% for i=1:size(data.trial,1)
%     for j=i:size(data.trial,1)
%         count=count+1; coeff_all(j,i) = coeff(count);
%         coeff_all(i,j) = coeff_all(j,i);
%     end
% end
% 
% coeff = mean(coeff_all,2);
% 
% [val,ref_trl] = max(coeff)
% % aligning with lv with xcorrelation
% cfg=[];
% cfg.data.trial(:,1,:) = auc; cfg.data.time = xax;% data.time all time extending beyond trial limit
% cfg.roi =[0-toi_len 1.5+toi_len]; cfg.method='topos'; cfg.classful=1; 
% cfg.shift = 200;
% cfg.ref_trl(1,:) = squeeze(data.trial(ref_trl,:,:));
% auc_aligned = lv_align(cfg);
% 
% xax_new = auc_aligned.time;
% auc_aligned=squeeze(auc_aligned.trial);
% ids = nearest(xax_new, 0):nearest(xax_new, 1.5);
% lv_pretty_errorbar(xax_new(ids), auc_aligned(:,ids), (auc_aligned(:,ids)*0)+0.5, 1);
