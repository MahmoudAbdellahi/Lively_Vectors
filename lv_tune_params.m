function varargout = lv_tune_params(varargin)
% many labels and values to be treated as defult values for the
% labels .. creates textarea for each.. and when the button is pressed it gets the new values or
% the default values if no changes ..

% EXAMPLE:
% [dimensions,blocks_perc,shift_perc,slide_trn] = lv_tune_params('do you want 1D classification?','0','percentage of trials per block','5',...
%             'shifting between blocks, percentage','1','are we sliding on train data?','1');

lbls = mod(1:numel(varargin),2);
vals_idx = find(lbls==0);
lbls_idx = find(lbls==1);


prompt = varargin(lbls_idx);
if any(ismember(prompt,'just use defaults'))
    params = varargin(vals_idx);
else
    dlgtitle = 'Set parameters';
    dims = [1 100];
    definput = varargin(vals_idx);
    params = inputdlg(prompt,dlgtitle,dims,definput);
end
% parameters = split(params)';
parameters=params;
for i=1:length(parameters)
    if ~isnan(str2double(parameters(i))), varargout{i,1} = str2double(parameters(i)); else % to autmatically convert what can be double to double
        varargout{i,1} = parameters(i);
    end
end




end



% old test for crewating dynamic positions
% function answer = lv_tune_params(varargin)
% % many labels and values to be treated as defult values for the
% % labels .. creates textarea for each.. and when the button is pressed it gets the new values or
% % the default values if no changes ..
%
% lbls = mod(1:numel(varargin),2);
% vals_idx = find(lbls==0);
% lbls_idx = find(lbls==1);
% % %
% % % position = [50 300 270 55];
% % % fig = uifigure;
% % %
% % % for i=1:length(vals_idx)
% % %     uilabel(fig,'Text',varargin{lbls_idx(i)},'Position', position );
% % %     % coordinates inside figure first two coordinates and the other two sizes of x and y of the textbox
% % %     position = position - [0 50 0 0];
% % %     txa(i) = uitextarea(fig,'HandleVisibility','on','Value',varargin{vals_idx(i)},'Position', position);
% % %     position = position - [0 50 0 0];
% % % end
% % %
% % %  f = figure
% % %  uibutton(fig, 'Text','Execute', 'ButtonPushedFcn', @(btn,event, f) btn_pressed(btn,txa, f));  % call back the save_fig function
%
% % if isfield(btn,'UserData'),
% %      btn.UserData
% %     updated_params = btn.UserData;
% % end
%
%
% prompt = varargin(lbls_idx);
% dlgtitle = 'Set parameters';
% dims = [1 35];
% definput = varargin(vals_idx);
% answer = inputdlg(prompt,dlgtitle,dims,definput);
%
%
%
%
%
% % updated_params = guidata(gcf);
%
% % fig = uifigure;
% % set(fig, 'Position',  [800, 800, 500, 400]) % the first two are the coordinates of the buttom corner and the second two are the size of the drawn figure in x,y
% %
% % assumed_path = strcat('D:\sul','''s code\Matt\sleep\erps\Organised\New exp\figures\mix');
% %
% %
% % txa_path = uitextarea(fig,'HandleVisibility','on','Value',assumed_path,'Position',[50 300 270 55]); % coordinates inside figure first two coordinates and the other two sizes of x and y of the textbox
% %
% % txa_folder = uitextarea(fig,'HandleVisibility','on','Value','Figure name','Position',[50 200 120 55]);
% %
% % btn = uibutton(fig, 'Text','Save', 'ButtonPushedFcn', @(btn,event) save_fig(btn,txa_path,txa_folder));  % call back the save_fig function
% % btn.Position = [50 100 150 50];
% end



% function  btn_pressed(btn,txa) % a call back function to save the figure when button is pressed
%
% length(txa)
% for i=1:length(txa)
%     updated_params{i,1} = txa(i).Value;
%
% end
% % btn.UserData = updated_params;
%
% % save data in gcf
%
% guidata(f, updated_params);
% % handles.updated_params = updated_params;
%
% % path = strcat(txa_path.Value, '\', txa_folder.Value);
% % saveas(gcf, strjoin([path '.emf']));
% % disp(strjoin(['Figure saved in: ' path '.emf']))
%
%
% end





