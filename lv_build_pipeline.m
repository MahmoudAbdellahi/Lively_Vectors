function result = lv_build_pipeline(cfg) 
% this function loads the cleaned data and starts with the steps of the
% pipeline as feature extraction and trials reduction and so on .. you can
% specify the sequence of execution in cfg.sequence and the path of data..
% and if there is no sequence then it's an automation and it will read the
% last record in automation that has no result and works on its sequence
% result is the final struct according to the last step of the sequence or
% the classification result if the seq. has classification at the end


if isfield(cfg,'sequence'), sequence = cfg.sequence; end
result = cfg.data;


func = sequence( mod(1:length(sequence),2) == 1);
values = sequence( mod(1:length(sequence),2) ~= 1);

method = cell(length(values),1); % initializing fields with empty and nans
window = nan(length(values),1);

% for lv_feature_extractor
id = find( ismember(func,'lv_feature_extractor') );
temp = string(split(values(id)));
method{id,1} = char(temp( mod(1:length(temp),2) == 1));% just incase we have the function executed many times 
window(id) = str2num(temp( mod(1:length(temp),2) ~= 1));



for i=1:length(func)
    cfg=[]; % all possible cfg fields are set here and if the functions doesn't have any of them it sould be nan or empty
    cfg.method = char(method(i)); 
    cfg.window = window(i);
    cfg.data = result;
    
    % functions with additional workarounds
    if strcmp(func(i),'lv_reduce_trials')
        good_trials = lv_reduce_trials(cfg);
        cfg2 = []; cfg2.trials = good_trials(1:no_good_trls); % when you use it you have to specify the no_good_trls so change it here
        result = ft_selectdata(cfg2, result);
    else
        func_handle = str2func(char(func(i))); % instead of eval to be able to debug code better
        result = func_handle(cfg); % the general case evaluate the function with the specified cfgs ... str2func
        %eval(['result = ' func(i) '(cfg);']); 
    end
    
    
    
end
 



end















