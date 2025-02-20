function mi = lv_mutual_information(cfg)
% takes data as observations x features and the labels of the two classes. Returns the mutual information
% for every feature
% Exhaustive on all observations for every feature .. uses entropy
% Generate data

% % example
% mean1 = [2, 2];
% cov = [1 0; 0 1];
% % mean2 = [-2, 6];
% mean2 = [-2, 2]; % having the same mean for feature 2 so that feature 1 is good
% rng(1)
% datax1= mvnrnd(mean1, cov, 100);
% datax2 = mvnrnd(mean2, cov, 100);  
% % Standardize data
% x = [datax1 ; datax2];
% y = [zeros(100, 1); ones(100, 1)]'; 
% 
% x = ((x - mean(x,1)) ./ std(x,[],1));
% 
% % Create a scatter plot
% figure;
% hold on;
% scatter(x(1:100, 1), x(1:100, 2), 'blue', 'o');
% scatter(x(101:end, 1), x(101:end, 2), 'red', 'x');
% xlabel('Feature 1');
% ylabel('Feature 2');
% legend('Class 1', 'Class 2');
% cfg.data = x; cfg.labels = y+1;


data = cfg.data;
labels = cfg.labels(:);

p1 = sum(labels==1)/length(labels);
p2 = sum(labels==2)/length(labels);

for feat=1:size(data,2)
    Ip = get_entropy(p1,p2);
    dat = data(:,feat);
    for obs=1:size(data,1)
        threshold = data(obs,feat);
        Dleft = sum(dat<=threshold); Dright = sum(dat>threshold);
        Dleft_c1 = sum(dat<=threshold & labels==1); Dleft_c2 = sum(dat<=threshold & labels==2);
        Ih_left = get_entropy(Dleft_c1/Dleft,Dleft_c2/Dleft);
        
        Dright_c1 = sum(dat>threshold & labels==1); Dright_c2 = sum(dat>threshold & labels==2);
        Ih_right = get_entropy(Dright_c1/Dright,Dright_c2/Dright);
        
        IG(obs,feat) = Ip - ((Dleft/length(labels))*Ih_left) - ((Dright/length(labels))*Ih_right);
    end
    lv_progress(feat,size(data,2),'Feature: ');
end


[mi,idx] = max(IG,[],1);

% hold on, 
% xline(data(idx(1),1));
% yline(data(idx(2),2));


end


function entropyVal = get_entropy(p1,p2)
entropyVal = -(p1*log2(p1+eps) + p2*log2(p2+eps));
end