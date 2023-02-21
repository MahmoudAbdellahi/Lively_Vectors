function [result] = lv_slider(cfg)
% returns: other_dim x samples x windows

% takes 3d or any dimensions and gets/compresses
% the sliding window: give it overlap and window size and it construct the
% windows in a higher dimension than the data (one dimension higher) and
% give it the higher dimension data and it will compress it .. a matrix of
% indices is used to know the original positions of windows and compressed
% data .. and if data.time is provided it will also be compressed along
% with the indices idx_marker ... the sliding and building of windows is
% done on the las dimenaion of the data the highest dim. so if 3d it will
% do windows on 3rd dim.

data = cfg.data;
shift = cfg.shift; % shifts in samples between windows (minimum is one sample shift because then we will get the full resolution)
windlen = cfg.windlen; % window length in samples

fprintf(['\n Sliding on data and getting windows, puts windows in higher dimension \n']);



sz = size(data.trial);

dim_of_interest = sz(end); % dimension of sliding
idx_marker = 1:dim_of_interest;

windows = (1:shift:dim_of_interest-windlen)';
if isempty(windows), result.higher_data=data.trial; return; end % if the window length matches data length then its the same as the data

for i=1:length(windows), temp(i,:) = [windows(i,1) : windows(i,1)+windlen]; end
windows = temp;
windows_as_vector = windows'; wind_sz=size(windows);

other_dim = repmat(':,',1,length(sz)-1); % to have the colon : for as many dimensions that we have except the last one (because it's of interest)
temp = eval(['data.trial( ' other_dim '  windows_as_vector(:));']); % it gets the indices of windows but after each other so we will need to reshape fot the windows to be 2d

other_dim = sz(1:end-1);
result.higher_data = reshape( temp, [other_dim wind_sz(2) wind_sz(1)] ); % other_dim x samples x windows

%         % to cells of windows
%         other_dim = num2str(sz(1:end-1)); % now actual sizes to be able to get them from data and convert to cell
%         other_dim = strrep(other_dim,'  ', ','); % replace spaces between numbers with comma(s) , for mat2cell to work
%         higher_data = eval(['mat2cell( temp, ' other_dim ', repmat(wind_sz(2),1,wind_sz(1)) );']); % grouping the actual windows in cells.. every window is a cell
%         % higher_data contains cells each cell has the data for every window .. so
%         % first cell has the multidimensional data of the first window and then
%         % second cell for second window and so on...

result.idx_marker = windows; % window_sample we can look at this one to know how the window was constructed and the idx of each element in the window

if isfield(data,'time'), result.time=data.time( round(median(windows,2)) ); end % to get the center index because that's the center of the window

fprintf('\n Done. Now the last dimension is replaced by: samples x windows. \n');


% example in rem wake: tested and working like the normal 100 smoothing
% cfg.data = im;
% cfg.shift = 1; % shifts in samples between windows (minimum is one sample shift because then we will get the full resolution)
% cfg.windlen = (win*200)/1000;
% result = lv_slider(cfg);
% size(result.higher_data)
% temp = squeeze(mean(result.higher_data,3));

end

