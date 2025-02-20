function [statsStruct,sval] = lv_mdim_clusterstats(cfg, result, parametric)
% mdim cluster permutation .. shuffle the sbj between cond and get the sample
% wise stats and sum the stats of contiguous points ... here we work on the
% difference between conditions so that if we want to shuffle a sbj
% between conditions we just flip the sign of the difference of the 'result'
% and compare to 0 and get the z-stats

% getting observed clusters
sample_alpha = cfg.clustercritval;
final_alpha = cfg.alpha;
permutations = cfg.n_permutations;

dims = repmat(':,',1,length(size(result{1,1}.perf))); dims=dims(1:end-1);

for i=1:length(result), eval(['cond_difference(i, ' dims ' ) = result{i,1}.perf;']); end
clear result;
sval = lv_mdim_pt_stats(cond_difference, parametric); 

connectivity = (length(size(cond_difference))-1)*2; % find the 4 connected points if 2d o.w. search in the mdim space
if length(size(cond_difference))==2, connectivity=1; end % incase of curve one connected point to the side
L = bwlabeln(sval>sample_alpha, connectivity); 
fprintf(['\n ' num2str(max(L(:))) ' cluster(s) found. \n']);


for i=1:max(L(:)), ob_cluster_sum(i) = nansum( reshape(sval.*(L==i),1,[]) ); end

decimal_val_lim = bi2de( ones(1,size(cond_difference,1)) ); % max decimal values if all sbj will be shuffled (all ones), then pick a random value then go binary to see who gets shuffled(1)
shuffle_info = de2bi( randi(decimal_val_lim, [1,permutations]) );

clusters_counters = zeros(1,length(ob_cluster_sum));% counter of the times the shuffled data is higher than the observed val.

% getting permutations
for i=1:permutations
    lv_progress(i,permutations,'performing cluster based permutation: ');
    shuff_cluster_sum=[];
    cond_difference_shuffled = cond_difference;
    eval(['cond_difference_shuffled(shuffle_info(i,:)==1,' dims ') = cond_difference(shuffle_info(i,:)==1,' dims ') .* -1;']); % flipping the sign to indicate shuffle
    sval_shuffled = lv_mdim_pt_stats(cond_difference_shuffled, parametric);
    
    L_shuff = bwlabeln(sval_shuffled>sample_alpha, connectivity);
    if max(L_shuff(:))>0
        for ii=1:max(L_shuff(:)), shuff_cluster_sum(ii) = nansum( reshape(sval_shuffled.*(L_shuff==ii),1,[]) ); end
        clusters_counters( max(shuff_cluster_sum) > ob_cluster_sum ) = clusters_counters( max(shuff_cluster_sum) > ob_cluster_sum )+1;
    end 

end 

statsStruct.p = clusters_counters ./ permutations;

statsStruct.mask = zeros(size(L));

idx = find(statsStruct.p<final_alpha);
for i=1:length(idx), statsStruct.mask(L==idx(i)) = 1; end % to be one despite the cluster

statsStruct.mask_with_cluster_numbers = statsStruct.mask .* L; % to show the actual cluster number

end