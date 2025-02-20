function transformed = lv_EEG_to_CNN(cfg)
% EEG topos to CNN activations, feature extraction..
% transforms data from 3d EEG using the training set to activation units of
% CNN after doing the convolution on the channels to reduce to 1d and this
% can then be used with RNN or other classifiers .. so it cares about the
% spatial proximity of channels ..
% returns the transformed trn and each transformed trn has a full tst
% because each trn time pt transforms the data with CNN weights..
trn = cfg.trn;
tst = cfg.tst;
transformed = [];

if ~isfield(trn,'label'), load lv_layout.mat; trn.label=lv_layout.label; end

cfg = []; % preparing the layout
cfg.layout = 'easycapM1.mat';
lay = ft_prepare_layout(cfg);
% cfg.rotate =90;
lay.label = lower(lay.label); labels = lower(trn.label);  % because the labels were sometimes capital and sometimes small
for i=1:length(labels)
    id = ismember(lay.label,labels(i)); temp_lay.label(i)=lay.label(id); temp_lay.pos(i,:)=lay.pos(id,:); temp_lay.height(i)=lay.height(id);
    temp_lay.width(i)=lay.width(id);
end
temp_lay.outline=lay.outline; temp_lay.mask=lay.mask;

% progressbar = ParforProgressbar(size(tst.trial,3),'title', 'Classification progress');
parfor trn_time=1:size(trn.trial,3)
    cfg=[]; % converting time pt to topo.
    cfg.data=trn;
    cfg.data.trial = cfg.data.trial(:,:,trn_time);
    cfg.method = 'eeg_to_topos_video';
    cfg.data.layout = temp_lay;
    peak_topo = lv_feature_extractor(cfg);
    % replacing nans with a baseline which is the mean of the image
%     baseline=nanmean(peak_topo.trial,2);
%     peak_topo.trial(:, isnan(mean(peak_topo.trial,1))) = repmat(baseline,1, sum(isnan(mean(peak_topo.trial,1))) );

    cfg=[]; % getting the model of the peak
    cfg.data.trial=peak_topo.trial;
    cfg.data.trialinfo = trn.trialinfo(:,1);
    cfg.transform = 0;
    cfg.method = 'CNN'; cfg.fs=200;
    cnn_model = lv_feature_extractor(cfg);
    % now we use the data to transform using the learned peak_topo cnn
    % so this should be again a loop on all topos of all trials and time
    % pts but now they will be reduced using cnn feature extraction

    % applying CNN weights to training set
    cfg=[];
    cfg.data_to_transform.trial = peak_topo.trial;
    cfg.transform = 1; cfg.net = cnn_model;
    cfg.method = 'CNN'; cfg.fs=200;
    temp = lv_feature_extractor(cfg);
    transformed_trial(:,:,trn_time) = temp.trial;


    % applying CNN weights to testing set
    for i=1:size(tst.trial,3)
        cfg=[];
        cfg.data=tst;
        cfg.data.trial = cfg.data.trial(:,:,i); % working on the current time pt
        cfg.method = 'eeg_to_topos_video';
        cfg.data.layout = temp_lay;
        topos_data = lv_feature_extractor(cfg);

        % replacing nans with a baseline which is the mean of the image
%         baseline=nanmean(topos_data.trial,2);
%         topos_data.trial(:, isnan(mean(topos_data.trial,1))) = repmat(baseline,1, sum(isnan(mean(topos_data.trial,1))) );

        cfg=[]; % transforming with the model of the peak
        cfg.data_to_transform.trial = topos_data.trial;
        cfg.transform = 1; cfg.net = cnn_model;
        cfg.method = 'CNN'; cfg.fs=200;
        temp = lv_feature_extractor(cfg);
        transformed_tst{1,trn_time}.trial(:,:,i) = temp.trial; % every trn_time has a transformed tst set
%         progressbar.increment();
    end

end
% delete(progressbar);
transformed.trn.trial = transformed_trial;
transformed.tst = transformed_tst;
