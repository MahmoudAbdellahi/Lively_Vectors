function [mdl] = lv_lda(cfg)
% LDA as in Matthias' paper (FDA) for binary classification and primal form
% Mdl = fitcdiscr(TRAIN,TRAIN_GROUP);
% [outclass,posterior] = predict(Mdl,TEST); err=1;
% takes mdl as is empty during training

% using mean difference and shirnkage regularisation .. using GED ... using
% the assumption of LDA and calculating the soft max of likelihood ..
%%%%
% cfg.use_GED=0;
mdl = cfg.mdl;
if strcmp(cfg.phase,'trn')
    if ~isfield(cfg,'use_GED')
        TRAIN = cfg.TRAIN;
        TRAIN_GROUP = cfg.TRAIN_GROUP;
        sz = size(TRAIN); n=sz(1);
        
        id1 = TRAIN_GROUP==1; id2 = TRAIN_GROUP==2;
        n1 = sum(id1); n2 = sum(id2);
        
        X = [ TRAIN(id1,:)-repmat(mean(TRAIN(id1,:),1),sum(id1),1); TRAIN(id2,:)-repmat(mean(TRAIN(id2,:),1),sum(id2),1)];
        X = X - repmat(mean(X,1),n,1);
        
        % Fisher discriminant analysis
        % Dual Ledoit-Wolf estimate of shrinkage parameter lambda that
        % controls the shape of the cov, 1 will be spherical and 0.5 is
        % spherical and so on
        S = (X'*X)/n; p1=0;
        for i=1:n, p1 = p1 + norm(X(i,:)' * X(i,:) - S,'fro')^2; end
        p2 = (n^2) * (trace(S^2)-(trace(S)^2)/sz(2));
        lambda = p1/p2; lambda = min(max(lambda,0),1); % to avoid overshooting and undershooting and always have lambda [0 1]
        
        Swithin = n1*(cov(TRAIN(id1,:),1)) + n2*(cov(TRAIN(id2,:),1));
        Swithin_hat = (1-lambda)*Swithin + lambda * trace(Swithin)/sz(2) * eye(sz(2));
        
        mdl.weights = (Swithin_hat^-1)*(mean(TRAIN(id1,:),1)' - mean(TRAIN(id2,:),1)');
        mdl.bias = -mdl.weights' * ((mean(TRAIN(id1,:),1)')+(mean(TRAIN(id2,:),1)'))/2; % the impact of bias is small
        % because based on that the sign will determine the class and since we took between class as mu1-mu2
        % now we assumed that class1 is positive but if the data is not centered the boundary will not be at zero
        % so we would need to shift it so that now the boundary at zero and the sign will reflect the class .. so we take
        % the central point of the classes (mu1+mu2)/2 and we use this as
        % our treshold .. so if the data is centered the effect of this
        % bias is small .. and we get the central point to do that because
        % it will be in the middle and should be reflecting zero so if it
        % is shifted that would be the bias..
        mdl.mus = [mean(TRAIN(id1,:),1); mean(TRAIN(id2,:),1)];
        mdl.Swithin_hat = Swithin_hat; mdl.sz=sz;
        mdl.cov1 = cov(TRAIN(id1,:),1); mdl.cov2 = cov(TRAIN(id2,:),1);
    else % GED gerenalised eigen decomposition
        TRAIN = cfg.TRAIN;
        TRAIN_GROUP = cfg.TRAIN_GROUP;
        sz = size(TRAIN); n=sz(1);
        
        id1 = TRAIN_GROUP==1; id2 = TRAIN_GROUP==2; 
         
        Sbetween =  cov( [mean(TRAIN(id1,:),1);mean(TRAIN(id2,:),1)]  ,1  ); 
        Swithin=  cov(TRAIN(id1,:),1) +  cov(TRAIN(id2,:),1); 
        if isfield(cfg,'tst_noise_for_within'), Swithin = cov(cfg.TEST,1); end % guess the true labels of tst should be revealed
        % to get the within cov correctly
        
        [W,val] = eig(Sbetween,Swithin); % generalized decomposition
        [~,id] = sort(diag(val),'descend');
        mdl.weights =  W(:,id(1));
        % to know the class direction
        % if class1 center is negative then we swap the weight sign because
        % this means that class1 is negative and classes directions are in
        % the other direction so that when we calculate mdl.outclass the result will be fine
        % we do that because we normally assume that class1 is positive at
        % mdl.outclass = 1*(Z>0)+2*(Z<0); so we do this to know the actual direction
        % this was not needed in the previous condition because the
        % difference of the mean was takes so the difference will be
        % positive when class1 cecnter is higheer than two so this handles the sign and gives direction to the result
        if mean(TRAIN(id1,:),1)*mdl.weights < 0, mdl.weights=-mdl.weights; end
        
        mdl.bias = -mdl.weights' * ((mean(TRAIN(id1,:),1)')+(mean(TRAIN(id2,:),1)'))/2; % the impact of bias is small
        
        mdl.mus = [mean(TRAIN(id1,:),1); mean(TRAIN(id2,:),1)];
        mdl.Swithin_hat = Swithin*n; mdl.sz=sz; % this *n will be removed when we divide by n later in prediction
        mdl.cov1 = cov(TRAIN(id1,:),1); mdl.cov2 = cov(TRAIN(id2,:),1);
    end
end

if strcmp(cfg.phase,'predict')
    TEST = cfg.TEST;
    Z = (TEST*mdl.weights) + mdl.bias;
    mdl.dval = Z;
    mdl.outclass = 1*(Z>0)+2*(Z<0);
    
    % posterior probability .. likelihood is from the ditribution 
    % then we softmax the likelihood .. the assumptions of likelihood are
    % that we have multivariate gaussian and the cov is the within cov and
    % is the same for both classes .. the difference between lda and bayes
    % classifier is that bayes classifier assumes independency of features
    % and get the pdf of every feature then multiply them to get the
    % likelihood for one class then softmax with other classes to get the
    % posterior .. but lda assumes multi-variate gaussian and same
    % covariance matrices for both classes
    likelihood1 = mvnpdf(cfg.TEST, mdl.mus(1,:), mdl.Swithin_hat/mdl.sz(1)); % dividing by
    %     the no. trials because the within cov was calculated for every trials and then summed 
    likelihood2 = mvnpdf(cfg.TEST, mdl.mus(2,:), mdl.Swithin_hat/mdl.sz(1)); % dividing by
    % softmax the likelihood to get the posterior
    mdl.posterior = [likelihood1./(likelihood1+likelihood2)  likelihood2./(likelihood1+likelihood2)];
    mdl.outclass = 1*(mdl.posterior(:,1)>mdl.posterior(:,2))+2*(mdl.posterior(:,1)<mdl.posterior(:,2));

   % showing that e^z isn't posterior for lda but for neural network is posterior 
%    post = zeros(size(TEST,1),2);
%    post(Z>0,1) = exp(Z(Z>0)); post(Z<0,2) = exp(Z(Z<0));
%    mdl.posterior = post./sum(post,2);

% %     % t-distribution lda not normal distribution
% %     df=mdl.sz(1)-1; % degree of freedom .. equals gaussian when df=inf.. and cauchy distribution when df=1 .. a parameter to optimise
% %     likelihood1 = mvtpdf(bsxfun(@minus,cfg.TEST,mdl.mus(1,:)),mdl.cov1,df); % we subtract the mean of the class because mvtpdf has the assumption that the t-dist. is at zero mean
% %     % and in equation 2.159 the mean was subtracted
% %     likelihood2 = mvtpdf(bsxfun(@minus,cfg.TEST,mdl.mus(2,:)),mdl.cov2,df);
% %     mdl.posterior = [likelihood1./(likelihood1+likelihood2)  likelihood2./(likelihood1+likelihood2)];
% %     mdl.outclass = 1*(mdl.posterior(:,1)>mdl.posterior(:,2))+2*(mdl.posterior(:,1)<mdl.posterior(:,2));
      
end






end

