

% from here: https://www.kaggle.com/datasets/broach/button-tone-sz?resource=download
% we can get the data from ERPdata.csv which contains the averages
% then we can use demographic.csv to see the labels and info .. then
% because the data are from each ppnt we could repeat the erp of each ppnt
% to have it like trials and then get the average of that and then use
% lv_erp both for averaging the trials and then for visualising the outcome
% on topos and doing the stats.. 


Data = readmatrix('D:\sul''s code\Matt\sleep\erps\Organised\New exp\lv_examples\ERP\data\1.csv\1.csv');