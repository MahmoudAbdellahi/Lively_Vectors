function cov_data_pooled=lv_pooled_cov(data)
% gets the pooled covariance of different trials (trials in cells)

data_centered = cellfun(@(x) (x-repmat(mean(x,2),1,size(x,2)) ), data,'Un',0);
cov_data = cellfun(@(x) ((x*x')./size(x,2)), data_centered,'Un',0);
cov_data_pooled=zeros(size(cov_data{1,1}));
for i=1:length(data), cov_data_pooled = cov_data_pooled + cov_data{1,i}; end

cov_data_pooled = cov_data_pooled./length(data);

end