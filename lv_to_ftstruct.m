function data_ft = lv_to_ftstruct(cfg)
% converts data astruct with different fields to a fieldtrip structure
% that's ready for fieldtrip functions, so it includes the necessary
% fields:
% .trial  .trialinfo  .sampleinfo  .time
% add to this functions more fields as necessary .. check they exist and
% then add them if not and maybe an option to have them or not ..
data = cfg.data;
timeax = data.time;

if ~isfield(data,'trialinfo'), data.trialinfo=ones(size(data.trial,1),1); end % assumed classless with all ones
if ~isfield(data,'sampleinfo') % assumed continuous sampleinfo with trials back to back
    temp = 1:length(timeax):length(timeax)*size(data.trial,1);
    data.sampleinfo = [temp' temp'+length(timeax)-1];
end
if ~isfield(data,'label') % assumed channel numbering just to have this field if labels are missing
    for i=1:size(data.trial,2)
        data.label{i,1} = ['channel' num2str(i)];
    end
end
if ~isfield(data,'dimord'), data.dimord='rpt_chan_time'; end
    
data_ft = data;

end



%% EXAMPLE
% load('sampleEEGdata.mat'); % good to use this to compare to his analyses
% data = []; 
% data.trial = permute(EEG.data,[3 1 2]);
% data.time = EEG.times/1000; 
% 
% cfg=[]; cfg.data=data;
% data_ft = lv_to_ftstruct(cfg);


