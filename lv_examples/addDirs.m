




% datadir = strcat('D:\sul','''s code\Matt\sleep\erps\Organised\New exp');

mri_append = ''; % for MRI add MRI_ and also check add Dirs because they now refer to different path o.w. keep empty
warning(['working on ' mri_append ' sleep data']);
if isempty(mri_append)
    rawdir = 'E:\jittered exp data before MRI/RawData/allSleep/data';
else
    % for MRI
    rawdir = 'E:\MRI_jittered_exp\allSleep\data';
end
 
% gender and age 
% before MRI
gender = [2 2 1 2 1 2 2 1 1 2 1 2 1 2 2 1 1 1 2 2 1  ] 
age = [21 19 28 18 22 20 22 22 24 20 18 19 20 18 19 18 21 21 19 20 22]
% after MRI
gender = [gender 1 2 1 2 2 2 1 2 2 2 1 2 1 2 2 2 1 1 1 1 1 1 2  1 2 2 1] 
age = [age 19 21 21 20 19 21 22 19 21 19 21 23 20 21 19 22 20 21 20 21 20 20 18 19 21 19 19]   

% ppnt code 18 low EHI score & no pre-sleep sequence learning .. male and
% age 21

% preprocdir = [rawdir '/cleaned/'];
% figdir= [datadir '/figures/']; 
% resultsdir= [datadir '/results/'];
% scoringdir = strcat([datadir '/Hypnos/']);

26 + 31
%% data from imagesc
% sw = getimage(gcf);
% 
% figure
% imagesc(sw)
% set(gca,'YDir','normal')
% 
% Sfig4c = sw;
% save Sfig4c Sfig4c;

%% 1d permutation
% 
% load('exp_auc.mat')
% load('adp_auc.mat')
% load('xax.mat')
% 
% cond1 = exp_auc;
% cond2 = adp_auc;
% timeax = xax;
% 
% idx = nearest(timeax, -0.04):nearest(timeax, 1.54);
% cond1 = cond1(:,idx);
% cond2 = cond2(:,idx); timeax = timeax(idx);
% mo_pretty_errorbar(timeax ,cond1,cond2, 99);
% 
% % cond1 and cond2 ... sbjxtime w el data gwaha erp maslan ..
% maskStat = nan(size(timeax));
% for i=1:length(timeax)
%     [P_wilcoxon,~,Wilcoxon_stats] = signrank(cond1(:,i), cond2(:,i), 'method' ,'approximate');
%     if P_wilcoxon<0.05
%         maskStat(i) = Wilcoxon_stats.zval;
%     end
% end
% % plot(timeax, maskStat)
% mask = ~isnan(maskStat);
% signiPts = bwconncomp(mask);
% signiIdx = signiPts.PixelIdxList;
% for i=1:length(signiIdx)
%     signiClusterStat(1,i) = sum( maskStat(signiIdx{1,i}));
% end
% 
% % dlwa2ty m3ak el clusters stats nbda2 el shuffling
% fprintf( '\n' );
% ft_progress('init', 'text',    'Please wait...')
% shuffles=1000; distrbutionShuffled=[]; 
% for i=1:shuffles
%     randshuff = randi([0 1],[1 size(cond1,1)]); 
%     temp1 = cond1(randshuff==1,:);
%     cond1(randshuff==1,:) = cond2(randshuff==1,:);
%     cond2(randshuff==1,:) = temp1;
%     
%     % same but for shuffled
%     maskStat = nan(size(timeax));
%     for ii=1:length(timeax)
%         [P_wilcoxon,~,Wilcoxon_stats] = signrank(cond1(:,ii), cond2(:,ii), 'method' ,'approximate');
%         if P_wilcoxon<0.05
%             maskStat(ii) = Wilcoxon_stats.zval;
%         end
%     end
%     
%     mask = ~isnan(maskStat);
%     signiPts = bwconncomp(mask);
%     signiIdx2 = signiPts.PixelIdxList;
%     for ii=1:length(signiIdx2)
%         signiClusterStatShuffled(1,ii) = sum( maskStat(signiIdx2{1,ii}));
%     end
%     
%     distrbutionShuffled(1,i) = max(signiClusterStatShuffled); % 3shan babos 3la positive clusters w bs 3shan kda bashof el max
%     ft_progress(i/shuffles); 
% end
% 
% for i=1:length(signiClusterStat)
%     Pval(1,i) = sum(distrbutionShuffled > signiClusterStat(1,i)) ./ shuffles    
% end
% 
% signiPart = nan(size(timeax));
% 
% clusterSingificanceThr = 0.05;
% signiPart(signiIdx{1,Pval<clusterSingificanceThr} ) = 0.6;
% hold on,
% plot(timeax, signiPart);


 

% % sleep align test
% tst_trls = sleepTheta.trial(:,7,:); % sleepTheta.trial: trlsxchxtime
% [Aligned_tst_trls] = align_erps(tst_trls, 2);
%
% m1=squeeze(Aligned_tst_trls );
% m2=squeeze(tst_trls );
% pretty_errorbar(sleepTheta.time, m1,m2,99 );
%
% % fixed small numbers align test
% s=[];
% ss1 = [zeros(1,10) 1:10 10:-1:1] +30;
% ss2 = [1:10 10:-1:1  zeros(1,10) ] +30;
% s(:,1,:) = [ ss1   ;  ss2  ];
%
% [Aligned_s] = align_erps(s, 2);
% m1=squeeze(Aligned_s );
% m2=squeeze(s );
% pretty_errorbar(1:size(m1,2), m1,m2,99 );
%
% figure,
% plot(1:length(ss1), [ ss1  ;  ss2])
%
% figure,
% plot(1:length(ss1), [ m1(1,:)  ;  m1(2,:)])
% % wake align test
% % 2o2af 3nd line 605 fe txt rem sleep wake classification
% sleepTFc1 = sleepTF.trial(clabel_hand_train==1,7,:);
% sleepTFc2  = sleepTF.trial(clabel_hand_train==2,7,:); %CP3
% m1=squeeze(sleepTFc1 );
% m2=squeeze(sleepTFc2 );
% pretty_errorbar(sleepTF.time, m1,m2,99 );
%
% hold on,
% sleepTFc1 = align_erps( sleepTF.trial(clabel_hand_train==1,7,:) , 2);
% sleepTFc2  = align_erps( sleepTF.trial(clabel_hand_train==2,7,:) , 2); %CP3
% m1=squeeze(sleepTFc1 );
% m2=squeeze(sleepTFc2 );
% pretty_errorbar(sleepTF.time, m1,m2,99 );



