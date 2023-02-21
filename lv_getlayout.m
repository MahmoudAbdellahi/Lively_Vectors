function data = lv_getlayout(data)
% gets the layout for the channels inside data.label and puts them in
% data.lv_layout

cfg = [];
cfg.layout = 'easycapM1.mat';
lay = ft_prepare_layout(cfg);
lay.label = lower(lay.label); labels = lower(data.label);  % because the labels were sometimes capital and sometimes small
for i=1:length(labels)
    id = ismember(lay.label,labels(i)); temp_lay.label(i)=lay.label(id); temp_lay.pos(i,:)=lay.pos(id,:); temp_lay.height(i)=lay.height(id);
    temp_lay.width(i)=lay.width(id);
end
temp_lay.outline=lay.outline; temp_lay.mask=lay.mask;
data.lv_layout=temp_lay;

end