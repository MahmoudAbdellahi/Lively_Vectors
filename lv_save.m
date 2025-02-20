function lv_save(path,data,varargin)
% uses h5 format to save the data and it seperates data.trial and saves it
% in h5 because h5 can compress matrices ..  and other_data is the other
% things in the struct like .time .label ... they are saved in another file
% with the same  name but .mat not h5
% TO USE:
% takes optional argument with the name of the matrix that has double
% because that's the field that has the data and is the biggest so when
% you put its name it will be saved in h5 and the other parts in other data
% if not specified then the normal save of matlab is used..

 
if ~isempty(varargin)
    array_name = varargin{1};
    % data.trial
    h5create([ path '.h5'] ,['/' array_name], eval(['size(data.' array_name ')']) );
    h5write([ path '.h5'],['/' array_name], eval(['data.' array_name]));
    
    % other data
    eval(['data.' array_name '=[];'])
    if isfield(data,'cfg'), data.cfg=[]; end
    other_data=data;
    save(path, 'other_data', '-v7.3');
else
    save(path, 'data', '-v7.3');
end
    
    
end


% %% manual trials 
% %% Quick save and load using h5 and dividing data into two segments
% % every sbj will be saved in two files one of them is h5 data set that
% % has data.trial and the other is a mat file with the rest of data fields
% % delete cfg. before hand if you can it saves a lot of memory and proc.time
% for nn=2:4
%   tic
%   load([preprocdir 'part' num2str(sbj(nn)) '_' type '_cleaned_N' num2str(sleep_stage)], 'cleaned_data'); % will be the new one from manual inspection
%   data = cleaned_data; clear cleaned_data;
%   toc   
%   
%   % data.trial
%   h5create([preprocdir '\final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_cleaned_h5_N' num2str(sleep_stage) '.h5'] ,'/trial',size(data.trial));
%   h5write([preprocdir '\final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_cleaned_h5_N' num2str(sleep_stage) '.h5'], '/trial', data.trial);
%   
%   % other data
%   data.trial=[]; data.cfg=[]; other_data=data;
%   save([preprocdir 'final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_cleaned_h5_N' num2str(sleep_stage)], 'other_data', '-v7.3'); 
%   
%   fprintf(['\n Done, subject: ' num2str(sbj(nn)) '\n']);
%   
%   % to read
% %   load([preprocdir 'final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_cleaned_h5_N' num2str(sleep_stage)], 'other_data');
% %   data=other_data; clean other_data;
% %   data.trial = h5read([preprocdir '\final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_cleaned_h5_N' num2str(sleep_stage) '.h5'],'/trial');
% end
%  