function aligned = lv_align_pipeline(cfg)
% 1. choose a reference without lagging the signals
% 2. align the signals to the reference signal.
% alignment is classful or classless depending on the signal and whether to
% align similar shapes only or similar and different (positive and
% negative correlations)
% choose: roi(trial duraion) and shift amount and give the uncut data  to
% example: cfg.roi =[0 1.1];  cfg.shift = 400; (400ms shift)
% the function..
if isfield(cfg,'classful'),warning('classful alignment!'); else warning('classless alignment!'); end
% classful means that postiive correlation is considered so the shape of
% signals should be the same but classless means that the shape could be
% opposite and the correlation could be negative and that's good for a
% signal alignment between classes because then we want the
% different/similar shapes to be aligned to each other because we have different classes

data = cfg.data;

coeff=[];
hold_full_signals = data;
data.trial = data.trial(:,:,nearest(data.time,cfg.roi(1)):nearest(data.time,cfg.roi(2))); % region of interest

for i=1:size(data.trial,1)
    for j=i:size(data.trial,1)
        c = reshape( normxcorr2( squeeze(data.trial(i,:,:)), squeeze(data.trial(j,:,:)) ),size(data.trial,2),[] );
        % only looking at the logical parts without edge artifacts of the 2D convolution, so we remove the padded convolution results
        c(:,[1:size(data.trial,3)-1 end-(size(data.trial,3)-1)+1:end])=[];
        c([1:size(data.trial,2)-1 end-(size(data.trial,2)-1)+1:end],:)=[];
        if isfield(cfg,'classful'), coeff = [coeff c]; else
            coeff = [coeff abs(c)]; % see the commented part above we just do it in a row for parfor to work
        end
        %        lag(i,j) = -lag(j,i);
    end
end
count=0; % refilling with the results from parfor
for i=1:size(data.trial,1)
    for j=i:size(data.trial,1)
        count=count+1; coeff_all(j,i) = coeff(count);
        coeff_all(i,j) = coeff_all(j,i);
    end
end

coeff = mean(coeff_all,2);

[val,ref_trl] = max(coeff)
% aligning with lv with xcorrelation
cfg.data = hold_full_signals; % data.time all time extending beyond trial limit
cfg.method='topos';
cfg.ref_trl = reshape(squeeze(data.trial(ref_trl,:,:)),size(data.trial,2),[]);

aligned = lv_align(cfg);


end