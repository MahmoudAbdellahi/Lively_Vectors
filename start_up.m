% % fieldtrip
% addpath('D:\codeHub\matlab\Toolboxes\fieldtrip-20190419') % fieldtrip-20190419 fieldtrip-20250107 fieldtrip-20250114
% ft_defaults
% 
% %classifier scripts 
% addpath(genpath([pathAppend 'sul''s', ' code\Matt'])) 
% cd([pathAppend 'sul''s', ' code\Matt\sleep\erps\Organised\New exp']) 
% 
% set(0,'defaultfigurecolor',[1 1 1])
% set(0,'defaultAxesFontSize',16)
% set(0,'DefaultAxesTitleFontWeight','normal');
% 
% set(0, 'defaultFigureUnits', 'normalized')
% 
% set(0, 'defaultFigurePosition', [0.0789 0.4931 0.6 0.3]) % cubric pc
% % set(0, 'defaultFigurePosition', [0.0789 0.4931 0.3 0.3]) % laptop
% 
% % set(0, 'defaultFigurePosition', [1.0789 0.4931 0.2184 0.2931])
% 
% % s = settings;
% % s.matlab.fonts.custom.editor.FontToUse.PersonalValue = "CustomFont";
% % s.matlab.fonts.custom.editor.Name.PersonalValue = "Courier New";
% % s.matlab.fonts.custom.editor.Size.PersonalValue = 10; % in point
% 
% % set(0, 'defaultFigurePosition', [0.0789 0.4931 0.3 0.3]) % laptop
% 
% 
% % addpath(genpath(['C:\Users\Mahmoud Eed\AppData\Roaming\Microsoft\Windows\Start Menu\Programs']))
 
% fieldtrip
addpath('D:\codeHub\matlab\Toolboxes\fieldtrip-20190419') % fieldtrip-20190419 fieldtrip-20250107 fieldtrip-20250114
ft_defaults
% going to lv path
addpath('D:\codeHub\matlab\lv');
cd('D:\codeHub\matlab\lv');
% figure position and size
set(0,'defaultfigurecolor',[1 1 1])
set(0,'defaultAxesFontSize',12)
set(0,'DefaultAxesTitleFontWeight','normal');
set(0, 'defaultFigureUnits', 'normalized')
set(0, 'defaultFigurePosition', [0.6 0.6 0.16 0.25]) % cubric pc [left, bottom, width, height]
% set(0, 'defaultFigurePosition', [0.0789 0.4931 0.3 0.3]) % laptop

plot(1:10, 1:10)
 %% for github update
% lv on git is the one in newexp
% We don't have to inti or Git remote everytime .. so after the first time add commit push 
% Git init
% Git add filenames
% Git commit -m “sss”
% Git remote add origin https://github.com/MahmoudAbdellahi/Lively_Vectors
% Git push -u origin master .. if needed could used --force instead of -u
% but if there are updates on the remote not on the local this would delete
% the files on the remote .. and there is a protection role on git delete
% it first to be able to do this
% getting the names of lv_ files to update them
% % names = dir('*lv*.m'); s=[];
% % for i=1:length(names), s = [s names(i).name ' ']; end

% Or from matlab go on the files then source control then add to git then commit then push
% or do branching and work on the same version in isolation from one another
  

%% for lv tutorials
% field trip 
% cd(['D:\sul''s', ' code\Toolboxes\fieldtrip-20190419'])
% ft_defaults
% 
% %classifier scripts 
% addpath(genpath(['E:\Lively Vectors'])) 
% cd(['E:\Lively Vectors']) 
% 
% 
% set(0,'defaultfigurecolor',[1 1 1])
% set(0,'defaultAxesFontSize',16)
% set(0,'DefaultAxesTitleFontWeight','normal');
% 
% set(0, 'defaultFigureUnits', 'normalized')
% set(0, 'defaultFigurePosition', [1.0789 0.4931 0.2184 0.2931])
