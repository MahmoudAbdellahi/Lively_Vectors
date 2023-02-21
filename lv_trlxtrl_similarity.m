function agg_similarity = lv_trlxtrl_similarity(cfg)
% takes cfg including trn and tst(trlxchxtime) and does RSA
% and other trlxtrl similarity like doing RSA classifier and also phase
% similarity between trials .. this depends on the function we use to call
% lv_dim_handle search Notes 8 for 'trl_trl similarity RSA maslan' for this
% because correlation will take dim1 as the thing to correlate and then all
% other dims if you feed a 3d matrix the result of correlation may be huge
% dimsData1xdimsData2 without the dim of interest (doi) so we need to
% divide the trials into blocks and apply this on block pairs and do all
% permutations of blocks

data1.trial = single(cfg.data1.trial);
data2.trial = single(cfg.data2.trial);


% trying to divide the trials into blocks but it didn't work at the end the
% result changed when we change the block size !
% % % the max array of memory at the moment we call the function memory but we
% % % will limit ourselves to half of this to do this operation
% % arrSz = memory;
% % max_single_elements = (arrSz.MaxPossibleArrayBytes)/4; % single takes 4 bytes
% % max_elements = max_single_elements*.2; % we talk 2/3 so we keep the memory free
% % 
% % % blocks division to fit in memory .. we divide trials, channels won't matter and let's keep time fixed
% % % assume that the data1.trial and data2.trial as x
% % % because x_time1 * trl_time2 = max_elements .. so
% % trls_per_block_temp = floor(  max_elements ./ (size(data1.trial,3)*size(data2.trial,1)*size(data2.trial,3))  );
% % if trls_per_block_temp<1, error('consider splitting the time into quadrants and repeat again!'); end % when the size
% % % of data2 is really big and no trials could be chosen from data1
% % 
% % trls_per_block = min(trls_per_block_temp, size(data1.trial,1));
% % 
% % if trls_per_block==trls_per_block_temp % many blocks
% %     it = 1:trls_per_block:size(data1.trial,1)-trls_per_block; it=it';
% %     it = [it it+trls_per_block-1];
% %     if it(end)<size(data1.trial,1), it = [it;it(end)+1 size(data1.trial,1)]; end % the remaining trials after last full block
% %     weight_last_block = (diff(it(end,:))+1) / trls_per_block; % because we take the mean of similarity from blocks we have to weight the last block differently so the
% %     % contribution to the mean is fair and based on the nummber of trials
% %     % within every block
% %     % weight_all_blocks = 1/size(it,1); weight_last_block = weight_last_block/size(it,1); % to verify this should add to one: (weight_all_blocks*(size(it,1)-1))+weight_last_block
% %     % if round((weight_all_blocks*(size(it,1)-1))+weight_last_block) ~=1, error('revise the weights of blocks because it doesn''t fairly add to one!'); end
% % else % all elements of dataset
% %     it = [1 size(data1.trial,1)];
% % end

cfg_hold = cfg;
it = (1:size(data1.trial,1))';
agg_similarity=zeros(size(data1.trial,3),size(data2.trial,3));
progressbar = ParforProgressbar(size(it,1),'title', 'RSA progress');
parfor i=1:size(it,1)
    cfg=cfg_hold;
    cfg.data2.trial = data2.trial;
    cfg.data1.trial = data1.trial(it(i),:,:);%(it(i,1):it(i,2),:,:);
    %     cfg.data2.trial = data2.trial(it(i,1):it(i,2),:,:); this is
    %     problematic because it will force you to make a second loop to take all block pairs
    % cfg.same_data contains the range if we want to apply on same data ..
    if isfield(cfg_hold,'same_data'), cfg.same_data = [it(i,1) it(i,2)]; end
     
    similarity = squeeze(nanmean(nanmean(lv_dim_handle(cfg),1), 3)); % mean of trials .. now it becomes timeData1xtimeData2
    
% %     if trls_per_block==trls_per_block_temp
% %         if i~=size(it,1)
% %             agg_similarity = agg_similarity + atanh(similarity);
% %         else
% %             agg_similarity = agg_similarity + (atanh(similarity).*weight_last_block);
% %         end
% %     else, agg_similarity = agg_similarity + atanh(similarity);
% %     end
    agg_similarity = agg_similarity + atanh(similarity);

    progressbar.increment();
    %lv_progress(i,size(it,1),'Progress: ');
end

delete(progressbar);
agg_similarity = agg_similarity./size(it,1);

end

function example
%%
% cfg.data1.trial = im.trial;
% cfg.data1.trial = im.trial;
rng(1)
cfg=[];
cfg.data1.trial =  (randn(100,4,10));
cfg.data2.trial =  (randn(20,4,30) );
cfg.dim_of_interest1 = 2; cfg.dim_of_interest2 = 2;
cfg.handle = 'corr'; cfg.arguments = {'type' , 'spearman'};
result  = lv_trlxtrl_similarity(cfg);
end

















