function comp = lv_component_analysis(cfg)
% performs different component based calculations on the data
% ICA, PCA, CSP .. these methods are different. ICA gets components and the
% dimensionality of the data is intact. PCA reduces dim. CSP depends on
% two classes and maximizes the distance between them (spatial firlter x
% will max. the var for one class and minimize it for the other) so when
% calling the function keep that in mind and give the proper parameters
% accordingly.

% this function better works with cell arrays of trial so we take the 3d
% data and convert it to cell arrays and call the functions

% cfg.data will include the data without the EOG channel and that what the code does here it removes the effect_chwhen you want to
% remove eye movements! because the point is to see the effect of eye
% movements on EEG channels.
% cfg.method includes the methods: % ex: 'runica'

% if you are doing this like ICA for eye mov. removal, then do it after
% interpolating the bad trials (after lv_clean_segmented) (after cleaning your data)

% Input: lv:example from runica implementation in ft which is actually from EEGlab
%    data     = input data (chans,frames*epochs).  , lv: so it puts trl_time together in one dimension
%               Note that if data consists of multiple discontinuous epochs,
%               each epoch should be separately baseline-zero'd using
%                  >> data = rmbase(data,frames,basevector); ,, lv: we
%                  should be fine because yes our data from diff. trials
%                  but we did deman in the beginning so they are on the
%                  same scale.

% in ft_componentanalysis, line 374 ft does the concatenation of trl_time
% into one dimension.. so the data becomes chx trl_time
% Example:
% cfg.data = data;
% cfg.method = 'runica';
% cfg.option = 'remove'; % can be none if we don't want to apply the transformation here and just want to get the components
% cfg.effect_ch = {'EOG1','EOG2'};
% comp = lv_component_analysis(cfg);


data = cfg.data;
method = cfg.method;

fprintf(['\n Performing component analysis with: ' method ' \n']);
% 3d to cell arrays
for i=1:size(data.trial,1), temp{1,i} = squeeze(data.trial(i,:,:)); end
data.trial=temp; clear temp;
if isfield(cfg,'data2') % this if we have 2datasets and want to pool them together most likely with PCA only
    for i=1:size(cfg.data2.trial,1), temp2{1,i} = squeeze(cfg.data2.trial(i,:,:)); end
    data.trial=[data.trial temp2]; clear temp2;
end


if isfield(cfg,'option') % for ICA
    option = cfg.option;
    effect_ch = cfg.effect_ch; % list of EOG ch. or motor in case of keep!
    % remove the effect channel so that we can see its effect on other channels
    cfg_tmp = []; cfg_tmp.channel = data.label(~ismember(data.label,effect_ch));
    data_to_run = ft_selectdata(cfg_tmp, data);

    comp_cfg = [];
    comp_cfg.method = method;
    [comp] = ft_componentanalysis(comp_cfg, data_to_run);

    figure('NumberTitle', 'off', 'Name', 'All components');
    st_cfg = [];
    st_cfg.component = 1:length(comp.label);       % specify the component(s) that should be plotted
    st_cfg.layout    = 'EEG1020.lay'; % specify the layout file that should be used for plotting
    st_cfg.comment   = 'no';
    ft_topoplotIC(st_cfg, comp)
else
    data_to_run = data;
    comp.trial={1};
end



% making sure no NaNs in data
tmp = cell2mat(comp.trial); tmp2 = cell2mat(data.trial);
if any(isnan(tmp(:))) || any(isnan(tmp2(:))), warning('lv: data contains NaNs interpolate that before component analysis.'); end

%% ICA
if strcmp(method,'runica')
    if strcmp(cfg.similarity,'mig_similarity') % ICA
        %% Correlate ICA with EOG
        % components in every trial are correlated with t he two EOG channels
        % (acrosstime in every trial) and then take the median of these
        % correlations and correct and see this vs permutation where the time
        % pts of EOG axtime are circularly shifted positively or neg. to make
        % shuffled data and then correlated with components the same way .. do
        % this many times and then compare the correlation(corrected with fiher to become normal) to permutation to be abel to get sensible
        % p-val. then apply two thresholds one to keep the components with 0.8 of cumulative correlation rho...
        % and the other to remove 0.3 of the no. of channels whichever will
        % give less components to be removed will be used .. the thing that I
        % don't like is that the correlation is across time which makes the
        % correlation always strong and significant vs permutation so it's
        % always up to the threshold to determine the components to reject and
        % almost always it will be 0.3*no. ch. so it will depend on the no. of
        % channels which shouldn't be the case because if we have 200 channels
        % why would we reject 60 components !! and also the 0.8 of cumulative
        % rhos will need lots of components so around 100 comp. to make 0.8 of
        % rhos and that's 100 component so the problem is when we have lots of
        % channels the behaviour may become strange ... also permutation of shuffling time pts is stringent on shuffled that�s why they
        % will not be shuffled and have the same instantaneous correlations in the eeg in general so they will not correlate with EOG and won�t be a control
        fprintf('    > computing correlation of ICA components and EOG: \n')
        st_datRep = data; nm_numPerm  = 200;
        nm_maxCumsum= 0.8;
        nm_maxPercCh= 0.3;
        tic
        vt_idEOG	= ismember(st_datRep.label,effect_ch);
        vt_corr     = cellfun(@(x,y) corr(x',y(vt_idEOG,:)'),...
            comp.trial,st_datRep.trial,'UniformOutput',false);

        mx_origCorr	= fn_fisherz(abs(cell2mat(vt_corr)));
        vt_medCorr  = median(mx_origCorr,2);

        mx_permCorr = nan(size(mx_origCorr,1),nm_numPerm);
        vt_perm     = randperm(numel(st_datRep.time),nm_numPerm);
        vt_idSig    = 1:numel(st_datRep.time);

        for kk = 1:nm_numPerm
            if mod(kk,2)
                nm_sign	= +1;
            else
                nm_sign	= -1;
            end

            vt_cId  = circshift(vt_idSig,nm_sign*vt_perm(kk));
            vt_corr	= cellfun(@(x,y) corr(x',y(vt_idEOG,vt_cId)'),...
                comp.trial,st_datRep.trial,'UniformOutput',false);

            mx_permCorr(:,kk) = median(fn_fisherz(abs(cell2mat(vt_corr))),2);
            lv_progress(kk, nm_numPerm, 'working on permutations: ')
        end

        vt_pValues	= repmat(vt_medCorr,1,nm_numPerm);
        vt_pValues  = sum(mx_permCorr >= vt_pValues,2)./nm_numPerm;

        vt_corrCh   = fdr_bh(vt_pValues);  % correct for multiple comparisons ... 

        [vt_medCorr,vt_cId]	= sort(vt_medCorr,'descend');

        vt_pValues 	= vt_pValues(vt_cId);
        vt_corrCh   = vt_corrCh(vt_cId);

        vt_cumsum   = cumsum(vt_medCorr)/sum(vt_medCorr);
        vt_cumsum   = vt_cumsum <= nm_maxCumsum;

        vt_pValues  = vt_pValues(vt_cumsum & vt_corrCh);
        vt_medCorr  = vt_medCorr(vt_cumsum & vt_corrCh);
        vt_cId      = vt_cId(vt_cumsum & vt_corrCh);

        nm_maxNumCh = round(nm_maxPercCh.*size(comp.trial{1,1},1));

        fprintf(' %i ICA components correlated, max %i selected - ',...
            numel(vt_cId),nm_maxNumCh)

        if numel(vt_cId) > nm_maxNumCh
            vt_cId	= vt_cId(1:nm_maxNumCh);
        end
        toc
        comp_to_reject = vt_cId;

        % visualizing the components
        % topo
        figure('NumberTitle', 'off', 'Name', 'Rejected with miguels correlation corrected and double thresholded 0.3channels and 0.8rho');
        st_cfg = [];
        st_cfg.component = comp_to_reject;       % specify the component(s) that should be plotted
        st_cfg.layout    = 'EEG1020.lay'; % specify the layout file that should be used for plotting
        st_cfg.comment   = 'no';
        ft_topoplotIC(st_cfg, comp)



    elseif strcmp(cfg.similarity,'lv_similarity') % ICA
        %% distance axtime and topo scores for every trial
        % for every trial: projecting the components on the data and then getting the
        % correlation between that projection topo and the proj. of effect ch.
        % .. and that's the topo way converted to score and then softmax and 10% thresholded,, the axtime way is the correlation
        % between each comp and the effect ch across time and the correlation
        % is converted into score as the diff of correlation because eye mov.
        % is negating each other on eog1 and eog2 then that is converted to
        % softmax and threshold of:(10% or more of rho) is used..

        effect_ch_id = find(ismember(data.label,effect_ch));
        non_effect_ch_id = find(~ismember(data.label,effect_ch));
        axtime_corr = cellfun(@(x,y) corr(x',y(effect_ch_id,:)','type','spearman'), comp.trial,data.trial,'Un',0);
        comp_proj = cellfun(@(x,y) (x * y(non_effect_ch_id,:)'), comp.trial,data.trial,'Un',0); % now each row has component topography
        eog_proj = cellfun(@(x,y) (x(effect_ch_id,:) * y(non_effect_ch_id,:)'), data.trial,data.trial,'Un',0); % now each row has eog topography
        % important: the actual topography is the vector of the component
        % itself transposed then multiplied by the cov matrix of the data
        % W'*cov check mike's min.9 in lec. 229 .. but I think no problem
        % since this is manual procedure and we look at the result before
        % rejecting anyway..
        topo_corr = cellfun(@(x,y) corr(x',y','type','spearman'), comp_proj, eog_proj,'Un',0); % correlating the topographies of each component topog. and the two eog topos

        trls_dis = median( cell2mat( cellfun(@(x) (     abs(diff(x,[],2))    ), axtime_corr ,'Un',0) ) ,2);
        softmax_dis = trls_dis ./ sum(trls_dis)
        to_reject = find(softmax_dis > 0.1)


        trls_dis = median( cell2mat( cellfun(@(x) (     abs(diff(x,[],2))    ), topo_corr ,'Un',0) ) ,2);
        softmax_dis = trls_dis ./ sum(trls_dis)
        to_reject2 = find(softmax_dis > 0.1)


        comp_to_reject = union(to_reject,to_reject2)
        % visualizing the components
        % topo
        figure('NumberTitle', 'off', 'Name', 'Rejected with lv axtime and projected topo. correlation with every trial');
        st_cfg = [];
        st_cfg.component = comp_to_reject;       % specify the component(s) that should be plotted
        st_cfg.layout    = 'EEG1020.lay'; % specify the layout file that should be used for plotting
        st_cfg.comment   = 'no';
        ft_topoplotIC(st_cfg, comp)


    end

    if isfield(data,'comp_to_reject') % if the component is already set and you know which to keep/reject
        comp_to_reject = data.comp_to_reject;
    end

    comp_to_reject_id = zeros(1,length(comp.label));
    comp_to_reject_id(comp_to_reject)=1;


    if strcmp(option,'keep') ==1
        % keep components correlating with channels (ex: motor) so low acc. is good
        cfg_comp = [];
        cfg_comp.component = find(comp_to_reject_id==0);
        data	= ft_rejectcomponent(cfg_comp, comp, data); % projection on data

    elseif  strcmp(option,'remove') ==1
        % remove components correlating with channels (ex: eye movement) so low acc. is bad
        cfg_comp = [];
        cfg_comp.component = find(comp_to_reject_id==1);
        data	= ft_rejectcomponent(cfg_comp, comp, data);
    end
    comp = data; % to return
    % visualizing the components
    % topo
    %         figure
    %         st_cfg = [];
    %         st_cfg.component = 1:numel(comp.label);       % specify the component(s) that should be plotted
    %         st_cfg.layout    = 'EEG1020.lay'; % specify the layout file that should be used for plotting
    %         st_cfg.comment   = 'no';
    %         ft_topoplotIC(st_cfg, comp)
    % axtime
    %         st_cfg = [];
    %         st_cfg.layout       = 'EEG1020.lay'; % specify the layout file that should be used for plotting
    %         st_cfg.viewmode     = 'component';
    %         st_cfg.allowoverlap	= 'yes';
    %         ft_databrowser(st_cfg,  comp)
end

%% PCA
if strcmp(method,'pca')
    % data.trial should have cells of trials each has channel_time .. if
    % your data is samples_features then put every sample vector of
    % features in a cell to simulate that time is 1 (so a sample cell will have
    % features_1).
    % here we get the eig vectors and values and return that so you can
    % pick the PCs that you want based on the variance and give them back
    % to this function to project the data using them
    if strcmp(cfg.step,'calculate')
        comp=[];
        data = cell2mat(data_to_run.trial)';
        data(isnan(sum(data,2)), :)=[]; % removing nans
        if cfg.centered==0
            data = data - repmat(mean(data,1),size(data,1),1); % removing the mean because the mean changes the location of the PC and data and we need
            % to be scale/value invariant and not change the orientation
        end
        [comp.eigVects,~,comp.eigVals] = pca(data); % gets eig vectors and values and they are sorted
        comp.eigVals = comp.eigVals./sum(comp.eigVals); % comp.eigVals is now explained variance so percentage of variance in the PC
        %according to variance descendingly, cols in eigVects
        % are the eigen vectors so to project data with pc1: eigVects(:,1)' * data' ... assuming data is samples_features
    end
    if strcmp(cfg.step,'transform')
        if cfg.centered==0
            temp = cell2mat(data.trial); temp = temp'; % removing the mean
            temp = (temp - repmat(nanmean(temp,1),size(temp,1),1))';
            temp = reshape(temp,size(data.trial{1,1},1),size(data.trial{1,1},2), length(data.trial));
            for i=1:size(temp,3), temp2{1,i} = squeeze(temp(:,:,i)); end
            data.trial = temp2;
        end
        comp.trial = cellfun(@(x) (cfg.eigVects(:,logical(cfg.chosen))' * x), data.trial,'Un',0); %cfg.chosen is logical that contains the eigVects to keep
        % this comp.trial is the projection of the data using the eigVects
        % chosen in cfg.chosen
        % go to myPCA for 2d visualization
    end
end

%% CCA
if strcmp(method,'cca')
    % canonical correlation analysis (CCA)
    % takes two datasets (observations x features) and aligns the features together .. by getting a
    % and b transformations .. the first column of a get the linear
    % transformation of the old features to get a new column that maximally
    % correlated with the first column of b transforming the old features
    % of the second dataset .. and that's how the first new feature is
    % calculated for both datasets and then it continues for all features
    % and it need the no. of rows to match so that we can perform the correlation .. 
    [a,b] = canoncorr(cfg.dataset1,cfg.dataset2); % datasets should be observations x features
    comp.dataset1 = cfg.dataset1*a; comp.dataset2 = cfg.dataset2*b;
end
%% SVD
if strcmp(method,'svd')
    % same input format as PCA
    dat = cell2mat(data_to_run.trial)';
    [U,S,V] = svd(dat); % so now we decomposed the rectangular data into: U S(singular values) and V
    % now we can reduce the noise by replacing some values in the diagonal of S
    % with low values with zeros .. and we can build the data matrix again
    % good for reducing noise and compressing information not dimensionality
    % and also building noise if needed

    % now getting the percentage of info in every value using S.. you should reduce info by zeroing some S
    info_perc = diag(S);
    info_perc = info_perc./sum(info_perc); % percentage of info contained in each S value
    info_cum = cumsum(info_perc);

    if strcmp(cfg.step,'calculate')
        % now you have the cumulative info percentage in the data and
        % should call the function again with the value of info to keep
        comp.info_cum = info_cum ;
    end
    if strcmp(cfg.step,'transform') || strcmp(cfg.step,'transform_noise')
        % go to myPCA for visualization
        chosen = zeros(length(info_cum),1);  chosen(1:nearest(info_cum, cfg.keep_info_threshold))=1;
        info_to_keep = chosen; % percentage of info to keep in the tranformed data


        if strcmp(step,'transform_noise'), info_to_keep = ~info_to_keep; end % to build noise by removing the keep_info_threshold from data
        % BEAWARE: that it's not actually noise in the common sense of
        % noise because it contains the most abstract representaion of
        % patterns in the data so still containing the pattern but very
        % abstract .. see cameraman example in myPCA

        new_S = diag( info_to_keep.*diag(S) );
        % reconstructing the data
        comp.trial = U*new_S*V';

        % to cell again to match the data format (cellarray of trials)
        comp.trial = mat2cell(comp.trial, size(comp.trial,1), repmat(length(data.time),1,length(data.trial)) );
    end
end



%% CSP
if strcmp(method,'csp')
    % the info of discrimination between classes should be already in the
    % variance before you use CSP .. what CSP does it help maximizing this
    % difference of variance between classes which leads to better
    % classification ...  it is normally used on band power as this is when
    % the discr. arises in the variance .. but you may use it with any data
    % as long as the useful info of discr. is the variance so instead of
    % channels will be features and time is 1

    % because CSP is different we calculate it and project the data .. CSP
    % is different because we cannot determine the no. of CSP filters to
    % keep and which ones until we project all filters on the data and see
    % the data and do something like CV on the projected data with
    % different no. filters each time to know which to keep and the
    % number.. so, call the function many times .. and choose from the
    % projected data the no. filters that you want in every iteration ..
    % normally CSP filters are good on both ends of the spectrum meaning
    % that the first and last filters and good and as you go towards the
    % middle they are less informative...
    way = 2;
    if strcmp(cfg.step,'calculate')
        %         S1 = cov( cell2mat( data.trial(data.trialinfo==1) )' );
        %         S2 = cov( cell2mat( data.trial(data.trialinfo==2) )' );
        S1 = pooled_cov(data.trial(data.trialinfo==1) ,'normal');
        S2 = pooled_cov(data.trial(data.trialinfo==2) ,'normal');
        comp=[];
        if way == 1
            [W,val] = eig(S1,S1+S2); % generalized decomposition (diagonalization) of data in S1 and S1+S2
        else
            P = whiten(S1 + S2, 1e-14);        % the average cov and get the whitening
            [B,~,~] = svd(P * S1 * P');  % whiten one of the classes and see its eigs
            W = B' * P;
            val = diag(1:size(W,1));
        end
        
        [comp.eigVals,comp.id]=sort(diag(val),'descend');
        [comp.eigVals2,comp.id2]=sort(diag(val),'ascend');
        comp.W = W;

    end
    if strcmp(cfg.step,'transform')
        comp=cfg.comp;
        W = comp.W(:,comp.id); 

        % project on data ..
        comp.trial = cellfun(@(x) (W'*x), data.trial,'Un',0);
    end
end



%% Source separation with generalised eigen decomposition (GED)
if strcmp(method,'source separation')
    % seperates the data from noise by finding the vector that maxmises the
    % distance between data and noise so that we supress the dimension with
    % the variance that is similar to the variance of noise and only keep
    % the dim with real variance so that when projected on data we will egt
    % the real spatial source of activation .. it isn't dim reduction but
    % source seperation so dim will remain the same but the real effect
    % will be localised .. but it depends because you can and should choose
    % the eigenVectors with the highest eigValues

    % if we have two datasets data and noise and the info is in the variance
    % then we can do a GED on them directly .. but if the info isn't in the vaariance
    % then we need to do between class/within
    % class maximisation to find the vector(s) that will maximally
    % separarte our classes .. GED or lda_beamformer
    noise = cfg.noise;
    if strcmp(cfg.step,'calculate')
        if cfg.info_in_variance == 1
            for i=1:size(noise.trial,1), temp{1,i} = squeeze(noise.trial(i,:,:)); end % noise to cell arrays
            noise.trial = temp;
            comp=[];
            S = cov( cell2mat( data.trial )' );
            R = cov( cell2mat( noise.trial )' );
            [W,val] = eig(S,R);
            [comp.eigVals,comp.id]=sort(diag(val),'descend');
            comp.W = W;
        else
            % maxi. between and mini. within class
            if strcmp(cfg.approach,'GED')
                for i=1:size(noise.trial,1), temp{1,i} = squeeze(noise.trial(i,:,:)); end % noise to cell arrays
                noise.trial = temp; clear temp;
                cov_within = pooled_cov(data.trial,'riemann') + pooled_cov(noise.trial,'riemann');

                temp=zeros(size(data.trial{1,1}));
                for i=1:length(data.trial), temp = temp+data.trial{1,i}; end % getting the mean of trials for between class cov
                data_mean = temp./length(data.trial); clear temp;

                cov_between = cov( [squeeze(mean(cfg.noise.trial,1))' ; data_mean'] );
                comp=[];
                [W,val] = eig(cov_between,cov_within);
                [comp.eigVals,comp.id]=sort(diag(val),'descend');
                comp.W = W;
            else
                % LDA beamformer should be similar to GED .. it calculates
                % and tranformers the data producing a source activity
                % since it compresses time dimension to get between class
                % effect we need to locate the effect in time beforehand
                % and only get the erp_pattern based on that .. then apply
                % on full lengths
                data_mean = squeeze(mean(cfg.data.trial,1));
                noise_mean = squeeze(mean(cfg.noise.trial,1));
                comp=[];
                erp_pattern = mean(data_mean-noise_mean, 2);% (:,150:end) % difference of classes erps is considered instead of cov_between .. and compressed in time
                %agg_data = [cfg.noise.trial ; cfg.data.trial];
                for i=1:size(noise.trial,1), temp{1,i} = squeeze(noise.trial(i,:,:)); end % noise to cell arrays
                noise.trial = temp; clear temp;
                agg_data=[];
                agg_data.trial = [data.trial noise.trial];
                % for within effect we get the whole aggregated
                % data because their covs will be summed inside beamformer
                [comp.W,comp.projected(:,1,:),comp.covariance] = LDAbeamformer(erp_pattern, agg_data  );

                temp = cell2mat(comp.projected.trial);
                temp = reshape(temp,size(comp.projected.trial{1,1},1),size(comp.projected.trial{1,1},2),[]);
                temp = permute(temp, [3 1 2]);
                comp.projected=temp; clear temp;

                comp.clabel = [ones(1,size(cfg.data.trial,1)) 2*ones(1,size(cfg.noise.trial,1))];
            end
        end
    end
    if strcmp(cfg.step,'transform')
        comp=cfg.comp;
        W = comp.W(:,comp.id);
        % project on data ..
        comp.trial = cellfun(@(x) (W'*x), data.trial,'Un',0);
    end
end





% back to 3d
if isfield(comp,'trial')
    temp = cell2mat(comp.trial);
    temp = reshape(temp,size(comp.trial{1,1},1),size(comp.trial{1,1},2),[]);
    temp = permute(temp, [3 1 2]);
    comp=temp; clear temp;
end




end







% helping function
function cov_data_pooled=pooled_cov(data, option)
% gets the pooled covariance of different trials (trials in cells)
% riemannian or normal averaging of covs
data_centered = cellfun(@(x) (x-repmat(mean(x,2),1,size(x,2)) ), data,'Un',0);
cov_data = cellfun(@(x) ((x*x')./size(x,2)), data_centered,'Un',0);
cov_data_pooled=zeros(size(cov_data{1,1}));

if strcmp(option,'normal')==1
    for i=1:length(data), cov_data_pooled = cov_data_pooled + cov_data{1,i}; end
    cov_data_pooled = cov_data_pooled./length(data);
elseif strcmp(option,'riemann')==1
    covariances = nan( size(data{1,1},1),size(data{1,1},1),length(data) );
    for i=1:length(data), covariances(:,:,i) = cov_data{1,i}; end
    cov_data_pooled = mean_covariances(covariances,'riemann');
end


end








