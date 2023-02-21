function lv_save_fig
% ui handle and save figure ... it will save the last selected figure because
% it's the one stored in get current figure (gcf).. so if we have many
% figures just mark on the one you want to save in order to save it .. you
% can browse and select the directory to save it to and it saves as .emf
% and .fig and also it remains with you so you can reuse the window many
% times
fig = uifigure;
set(fig, 'Position',  [800, 800, 500, 400]) % the first two are the coordinates of the buttom corner and the second two are the size of the drawn figure in x,y

txa_folder = uitextarea(fig,'HandleVisibility','on','Value','Figure name','Position',[100 220 270 55]);% coordinates inside figure first two coordinates and the other two sizes of x and y of the textbox

btn = uibutton(fig, 'Text','Save as', 'ButtonPushedFcn', @(btn,event) save_fig(btn,txa_folder));  % call back the save_fig function
btn.Position = [100 100 150 50];
end



function   save_fig(btn, txa_folder) % a call back function to save the figure when button is pressed

[pathname] = uigetdir(pwd);

path = strcat(pathname, '\', txa_folder.Value);
saveas(gcf, strjoin([path '.emf']));
saveas(gcf, strjoin([path '.fig']));
disp(strjoin(['Figure saved in: ' path '.emf and .fig']))
end

