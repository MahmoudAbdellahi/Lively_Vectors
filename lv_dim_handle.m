function result = lv_dim_handle(cfg)
% takes dataset or two datasets and performs the function specified in
% cfg.handle on a spcific dim of interest and optional handle inputs..
% max dim of data is 4d


if isfield(cfg,'arguments'), arguments = cfg.arguments; end % optional function handle arguments that will appended
% as string during the handle call

if isfield(cfg,'data')
    dim_oi = cfg.dim_of_interest;
    data = cfg.data; % ones data set
    % to dim of interest in dim1 ..
    dims = 1:length(size(data.trial)); sz=size(data.trial); sz(dim_oi)=[];
    dims(1) = dims(dim_oi); dims(dim_oi)=1; % to swap the doi with dim1
    data.trial = permute(data.trial, dims); data.trial = reshape(data.trial,size(cfg.data.trial,dim_oi),prod(sz));
    % apply function handle
    if ~isfield(cfg,'arguments')
        eval(['result = ' cfg.handle '(data.trial);']); % handle e.g., @mean
    else
        eval(['result = ' cfg.handle '(data.trial,' arguments ');']); % handle e.g., @mean
    end
    % reshaping back
    result = reshape(result,[1 sz]); % 1 because doi is compressed after the handle call
    result = permute(result, dims);
    
else
    % this else is for similarity so expected output from the handle is dimsData1xdimsData2 .. dims here is the mutiplication of all other dims than doi
    data1 = cfg.data1; data2 = cfg.data2;
    dim_oi1 = cfg.dim_of_interest1; dim_oi2 = cfg.dim_of_interest2;
    % to dim of interest in dim1 ..
    dims = 1:length(size(data1.trial)); sz1=size(data1.trial); sz1(dim_oi1)=[];
    dims(1) = dims(dim_oi1); dims(dim_oi1)=1; % to swap the doi with dim1
    data1.trial = permute(data1.trial, dims); data1.trial = reshape(data1.trial,size(cfg.data1.trial,dim_oi1),prod(sz1));
    % to dim of interest in dim1 ..
    dims = 1:length(size(data2.trial)); sz2=size(data2.trial); sz2(dim_oi2)=[];
    dims(1) = dims(dim_oi2); dims(dim_oi2)=1; % to swap the doi with dim1
    data2.trial = permute(data2.trial, dims); data2.trial = reshape(data2.trial,size(cfg.data2.trial,dim_oi2),prod(sz2));
    
    % apply function handle
    if ~isfield(cfg,'arguments')
        eval(['result = ' cfg.handle '(data1.trial,data2.trial);']);
    else
        eval(['result = ' cfg.handle '(data1.trial,data2.trial,''' arguments{1} ''', ''' arguments{2} ''');']); % correlation name value pair
    end
    if isfield(cfg,'same_data')
        idx_block = 1:size(cfg.data1.trial,1); idx_data2=cfg.same_data(1):cfg.same_data(2); % the symmetric part because we take only the blovk of data1 and all trials of data2
        for j=1:length(idx_block), result(idx_block(j),idx_data2(j))=nan; end
    end% if same data then we should remove the repeating parts because trials will correlate with themselves
    % reshaping back
    result = reshape(result(:), [sz1 sz2]); % because result is dimsData1xdimsData2 dims are dims other than the doi
end



end








%% example run with evaluate selection
function example
%% data
x = randn(100,4,200);
cfg=[];
cfg.data.trial = x;
cfg.dim_of_interest = 2;
cfg.handle = 'mean';
result = lv_dim_handle(cfg)
% mean(x)
%% data1 and data2
rng(1)
cfg=[];
cfg.data1.trial =  (randn(100,4,10));
cfg.data2.trial =  (randn(20,4,30) );
cfg.dim_of_interest1 = 2; cfg.dim_of_interest2 = 2;
cfg.handle = 'corr'; cfg.arguments = {'type' , 'spearman'};
result  = lv_dim_handle(cfg);

% run manually and compare:
mx=[];
for i=1:size(cfg.data1.trial,1)
    for j=1:size(cfg.data2.trial,1)
        for t1=1:size(cfg.data1.trial,3)
            for t2=1:size(cfg.data2.trial,3)
                mx(i,t1,j,t2) = corr( squeeze(cfg.data1.trial(i,:,t1))',squeeze(cfg.data2.trial(j,:,t2))','type' , 'spearman' );
            end
        end
    end
end
hh = (result(:)~=mx(:));
wrong = find(hh==1); wrong(1) % fine no difference and that is tested and it's fine
end



