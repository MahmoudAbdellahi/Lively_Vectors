%% running a pipeline .. the new one in lv_pipeline
sequence = {'lv_feature_extractor','erp 100'}; % to run many sequences back to back put new seq. in rows so it will be run x sequence

if ~exist('sequence','var') % if the sequence not provided then load from excel sheet
    table_specs = readtable([pwd '\lv_automation_runs.xlsx']);
    titles = table_specs.Properties.VariableNames;
    variables = table_specs.Variables;
    run_id=[];
    for i=1:size(variables,1), if strcmp(char(variables(i,end)),'done')==0, run_id=[run_id ; i]; end, end
    sequence = variables(run_id,2:end);
    fprintf(['\n lv: Applying pipeline sequence(s): \n']); disp(sequence);
    record=repmat({' '},size(variables,1),1); record(run_id,1)={'done'};
    table_specs = [table_specs record];
    writetable(table_specs, [pwd '\lv_automation_runs.xlsx']); % will put done here better than checking again at the end of execution
end

pathname = ''; % path to save the figures of results
for run=1:size(sequence,1)
    for cond=1:conditions % example: exp, adp
        for nn=1:numel(sbj)
            trn = lv_load([preprocdir '\final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_cleaned_h5_N' num2str(sleep_stage)], 'trial');
            tst = lv_load([preprocdir '\final_cleaned_after_inspection\part' num2str(sbj(nn)) '_' type '_cleaned_h5_N' num2str(sleep_stage)], 'trial');
            
            cfg=[];
            cfg.data=trn.data;
            cfg.sequence = sequence;
            result_trn = lv_build_pipeline(cfg); % for trn
            
            cfg.data=tst.data;
            result_tst = lv_build_pipeline(cfg); % for tst
            
            % classification
            result_trn.trial = zscore(result_trn.trial,[],1); result_tst.trial = zscore(result_tst.trial,[],1); % zscoring
            cfg=[];
            cfg.method = 'timextime';
            cfg.classifier_type = 'lda';
            cfg.trn = result_trn; %has labels in trialinfo
            cfg.tst = result_tst; cfg.folds=nan;
            classification_result{nn,cond} = lv_classify(cfg);
            
        end
    end
    result_exp = classification_result(:,1);
    
    if vschance==1, result_chance = classification_result(:,1); else, result_chance = classification_result(:,2); end % chance or the control condition
    
    dat=[];
    for i=1:size(result_exp,1), cond1(i,:,:)=result_exp{i,1}.perf; cond2(i,:,:)=result_chance{i,1}.perf; end
    id = mod(1: size(cond1,1)+size(cond2,1), 2);
    if vschance==1, dat.trial(id==0,1,:,:) = (cond2.*0)+0.50; else, dat.trial(id==0,1,:,:) = cond2; end
    dat.trial(id==1,1,:,:) = cond1;
    dat.time=tst.data.time; dat.freq=trn.data.time; dat.label={'z-stat'};
    [ TF_struct ] = lv_tf(dat,1,1); %data, do_stats, do_plot
    
    
    % the gcf has the result now save it for every run,, la de figure 100 w
    % kda el ids w b3d kda close el figures ba3d ma tesave ..
    if ~isempty(run_id), run_name = ['excel_' num2str(rund_id(run))]; else, run_name=num2str(run); end % if read from excel use the id to name the figures
    path = strcat(pathname, '\', ['run ' run_name]);
    saveas(gcf, strjoin([path '.emf']));
    saveas(gcf, strjoin([path '.fig']));
    
    fprintf(['\n Finished run: ' run_name ' \n']);
end


%% axspace classifier .. careful before using this make sure the channels order in lay match the data
cfg=[];
cfg.method = 'axtime';
cfg.classifier_type = 'lda';
cfg.do_parallel = 1;
sleep.trial = permute(sleep.trial, [1 3 2]); % permute space and time
im.trial = permute(im.trial, [1 3 2]);

cfg.trn = sleep; cfg.trn.trialinfo=clabel_hand_test;
cfg.tst = im; cfg.folds=nan; cfg.tst.trialinfo=clabel_hand_train;
[ sxsAuc{nn,1}.perf ] = lv_classify(cfg); % 1d curve in space
% % % % % % da  b2a bara el loop na5od el mean w nshof 3la el topo el donia 3amla ezay aw nethreshold kol topo gwa el loop le kol sbj !
cfg = [];
cfg.layout = 'MO_layout.lay'; % looks like easycap as well !
anne_layout = ft_prepare_layout(cfg);
% cfg = [];
% cfg.layout = lay;   % this is the layout structure that you created with ft_prepare_layout
% ft_layoutplot(cfg);
ss = zeros(1,15); ss(4:6)=1;
ft_plot_topo(anne_layout.pos(:,1), anne_layout.pos(:,2), ss, 'mask', anne_layout.mask,'outline' ,anne_layout.outline); % 2rsem 3leha el zvalues!

%%  cleaning sleep with lv_clean_segmented and interpolation in space .. careful before using this make sure the channels order in lay match the data
cfg = [];
cfg.layout = 'MO_layout.lay'; % looks like easycap as well !
anne_layout = ft_prepare_layout(cfg);
sleep.lv_layout = anne_layout;
sleep.aoi = sleep.label;
sleep.trialinfo=clabel_hand_test;
stats_window = [sleep.time(1) sleep.time(end)];

cleaned_data = lv_clean_segmented(sleep, [], stats_window, sbj(nn));


%% feature explorer
% feature explorer feels different but actually we should just look at the
% dimensions first to find that we can use the TF plotting to show the
% feature explorer because TF is sbj_ch_freq_time and explorer is
% sbj_ch_trl_time so we can use it in this way .. will start by sorting the
% trials based on the measure and then call the function

% here we want the conditions below each other in the same plot so I won't
% put them in sbj I will repeat sbj twice ..

% example with variance high to low
data = im;
data.trialinfo = clabel_hand_train;

[trl_measure,trl_idx] = sort(median( var(data.trial,0,3), 2), 'descend'); % IMPORTANT: will appear ascend when plotted because imagesc does: ynormal
% .. also the classes will be C2 then C1, you can use:
% set(gca,'YDir','reverse'); to flip it and it will make sense
data.trial = data.trial(trl_idx, :,:); data.trialinfo = data.trialinfo(trl_idx);
data.trial = data.trial([find(data.trialinfo==1)';find(data.trialinfo==2)'],:,:); trl_measure=trl_measure([find(data.trialinfo==1)';find(data.trialinfo==2)']);
data.trialinfo = data.trialinfo([find(data.trialinfo==1) find(data.trialinfo==2)]);

% appending the measrue to the data
sz = size(data.trial);
fs=200;
data.trial = cat(3, data.trial, repmat(trl_measure, 1,sz(2), round(sz(3)/4)) );
data.time = [data.time data.time(end)+1/fs:1/fs: data.time(end)+((round(sz(3)/4))/fs) ];


cfg = [];
cfg.layout = 'easycapM1.mat';
lay = ft_prepare_layout(cfg);
select = ismember(lower(lay.label),  lower(data.label));
lay.pos = lay.pos(select,:); lay.height = lay.height(select); lay.label = lay.label(select); lay.width = lay.width(select);
% reordering data to match the layout ordering
newIdx=[];
for i=1:length(data.label), newIdx = [newIdx ; find( ismember( data.label , lay.label{i} ) )]; end
data.label = data.label(newIdx);
data.trial = data.trial(:,newIdx,:);

% just plotting the layout
cfg = [];
cfg.layout = lay;   % this is the layout structure that you created with ft_prepare_layout
ft_layoutplot(cfg);
set(gcf,'Name',['layout sbj:' num2str(sbj)],'NumberTitle','off');


% inside nn main loop
close all
cond1(1,:,:,:)  =  permute(data.trial,[2 1 3]); % sbj_ch_trl_time, sort .trial according to the measure and put them cond1 and below it cond2
dat=[];
id = mod(1: size(cond1,1)+size(cond1,1), 2);
dat.trial(id==0,:,:,:) = cond1.*0;
dat.trial(id==1,:,:,:) = cond1;
dat.time=data.time; dat.freq=1:size(data.trial,1); dat.label=data.label;
dat.lv_layout = lay;
[ TF_struct ] = lv_tf(dat,0,1); %data, do_stats, do_plot

figure(1000)
set(gca,'YDir','reverse');
caxis([-50 50])


%% lv n100 or phenomenon one class classifier
% example is in PHASE_analysis_hilbert
cfg=[];
cfg.method = 'n100';
cfg.data = data_agg; % array of structs containing 3d with only the channel(s) containing the phenomenon, data_agg{1,nn} = im.trial(:, ismember(lower(im.label),lower('cz')), :);
cfg.time = im.time;
cfg.routine = 'training';
[ good_trials ] = lv_reduce_trials(cfg);

% assumning testing set using one sbj
% WARNING: you have to use [0.070 0.120] with the same no of samples when
% you test so make sure to do that with new datasets because that's how the
% classifier was trained .. 11 time pts in [0.070 0.120]sec.
cfg=[];
cfg.method = 'n100';
cfg.data = data_agg{1};
cfg.time = im.time;
cfg.routine = 'testing';
[ good_trials ] = lv_reduce_trials(cfg);

%% python in matlab
% put the path to python so that matlab can load it .. it only need to load
% it once .. if you want to change the version use: pyversion(version)..
% then import the module that you want and use it

version = ('C:\Users\Mahmoud Eed\AppData\Local\Programs\Python\Python36\python.exe');
pyversion % use : pyversion(version) .. if you want to change to another version

list = py.list({'This','is a','list'})


mod = py.importlib.import_module('core');
mod = py.importlib.import_module('linalg');

x = py.numpy.linspace(0,10,101);


% install numpy
python setup.py install --install-lib D:\numpy\numpy-master

pysys = py.sys.path;
pysys.append('D:\numpy\numpy-master\Packages\numpy-1.20.0.dev0+unknown-py3.7-win-amd64.egg')



module = py.importlib.import_module('io');



pysys = py.sys.path;
pysys.append('D:\numpy');

module = py.importlib.import_module('D:\numpy\numpy-master\numpy\matrixlib\tests\test_numeric');


%
% Specify Python Executable Library.
% pcPythonExe = 'C:\Users\Mahmoud Eed\AppData\Local\Programs\Python\Python36\python.exe';
[ver, exec, loaded]	= pyversion;
% Ensure python-matlab integration code is on matlab path.
matlabroot = 'D:\matlab 2018a installed';
pyFolder = fullfile(matlabroot, 'toolbox', 'matlab', 'external', 'interfaces', 'python');
addpath(pyFolder);
insert(py.sys.path, int64(0), pyFolder);
% Add relevant libaries' directory to python path.
otherpyFolder = 'C:\Users\dmattioli\.PyCharm2019.1\system\python_stubs\278535617';
insert(py.sys.path, int64(0), otherpyFolder)


pyversion

% import numpy

test0	= py.importlib.import_module('numpy') % Succeeds
test1	= py.importlib.import_module('matplotllib') % Succeeds
test2	= py.importlib.import_module('tensorflow') % Succeeds
test3	= py.importlib.import_module('keras') % Fails



modpath = '/usr/lib/python2.7/dist-packages/fabio';

hh = 'C:\Users\Mahmoud Eed\AppData\Local\Programs\Python\Python36\Lib\site-packages\numpy';
% follow the approach from the matlab documentation center
P = py.sys.path;
if count(P,hh) == 0
    insert(P,int32(0),hh);
end
np = py.importlib.import_module('test');
np.array(10)

% nas kteer 2alo feha mashakel w even lw kara2 numpy msh madmon eno
% yshta3'al w kaman aslan fe libraries tanya tensorflow w keras w kda fa da
% kda msh madmon neha2e ...

np = py.importlib.import_module('test');

%% da sha3'al bs start from anaconda
% to setup package ro7 hena lel pip mn cmd : cd C:\Users\Mahmoud Eed\AppData\Local\Programs\Python\Python36\Scripts
% b3d kda pip install xxxx
np = py.importlib.import_module('numpy');
x = np.array([10,1])

x = np.arange(15)

reshape(x,[3 5 ])

keras = py.importlib.import_module('keras');
tf = py.importlib.import_module('tensorflow');

% E:\test tensorflow\mypython\Scripts\activate

%% another way just import the module that has the low lvl python function
func2 = py.importlib.import_module('my_hello_arithm');
hh = func2.calcu(500);

% my_hello_arithm.py is a python file with: .. so you can write the low lvl
% python function in python and call it in matlab to be able to transfer
% variables
% import numpy as np
% def calcu(c):
%  a = 10 #np.arange(15).reshape(3, 5);
%  z=a+c
%
%  return z



%% lv_component_analysis
%% ICA
% from fieldtrip prepare the data with EOG to use ICA on it
cfg = [];
cfg.dataset            = 'ArtifactMEG.ds';
cfg.trialdef.eventtype = 'trial';
cfg = ft_definetrial(cfg);
cfg.channel            = 'MEG';
cfg.continuous         = 'yes';
data = ft_preprocessing(cfg);
% downsample the data to speed up the next step
cfg = [];
cfg.resamplefs = 300;
cfg.detrend    = 'no';
data = ft_resampledata(cfg, data);

% cfg        = [];
% cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
% comp = ft_componentanalysis(cfg, data);
% figure
% cfg = [];
% cfg.component = 1:20;       % specify the component(s) that should be plotted
% cfg.layout    = 'CTF151.lay'; % specify the layout file that should be used for plotting
% cfg.comment   = 'no';
% ft_topoplotIC(cfg, comp)

cfg=[];
cfg.data = dat;
cfg.method = 'runica';
cfg.option = 'remove'; % can be none if we don't want to apply the transformation here and just want to get the components
cfg.effect_ch = {'EOG1','EOG2'};
cfg.similarity = 'mig_similarity';
% cfg.similarity = 'lv_similarity';
comp = lv_component_analysis(cfg);


%% PCA
% claculate
cfg=[];
cfg.data = im;
cfg.method = 'pca';
cfg.step = 'claculate';
comp = lv_component_analysis(cfg);

% transform
cfg.eigVects = comp.eigVects;
cfg.chosen = cumsum(comp.eigVals)<=0.8;
cfg.step = 'transform';
tranformed_data = lv_component_analysis(cfg);

% for example on PCA on trials check REM_ERP_txt_sleepWake



%% testing lv_phenomenon
cfg = [];
cfg.channel = 'Cz';
[data] = ft_selectdata(cfg, dat);

% % cfg = [];
% % cfg.channel = 'Cz';
% % [phenom_dat] = ft_selectdata(cfg, phenom_dat);
% from trials to contiunous
% dat_temp = data;
% dat_temp.dimord ='chan_time';
% dat_temp = rmfield(dat_temp,[{'trialinfo'} {'sampleinfo'}]);
% dat_temp.trial = mat2cell( cell2mat(dat_temp.trial), 1,length(cell2mat(dat_temp.trial))) ;
% dat_temp.time = mat2cell(0:1/dat_temp.fsample:(length(dat_temp.trial{1})-1)/dat_temp.fsample ,1,length(dat_temp.trial{1}) ); dat_temp.cfg=[];
% data = dat_temp;

cfg=[];
cfg.phenomenon = 'SO';
cfg.fun_seq = [{'band'} {'duration_zcrossings'} {'spec_thresholded'}];
cfg.data.fs = 200;
cfg.data = data;

[phenom] = lv_detect_phenomenon(cfg)

%spindles
cfg=[];
cfg.phenomenon = 'spindle';
cfg.fun_seq = [{'band'}  {'spec_thresholded'}];
cfg.data.fs = 200;
cfg.data = data;

[phenom] = lv_detect_phenomenon(cfg);



% visualizing the phenomenon from the unfiltered raw signal
data_visual = data;
all_idx = cell2mat(phenom.measures(:,1)'); % all phenom. indices
temp = zeros(size(data_visual.trial{1}));
temp(all_idx) = data_visual.trial{1}(all_idx);
data_visual.trial = mat2cell(temp, 1, length(data_visual.trial{1}));
artf=ft_databrowser([],data_visual);


% visualizing locking on SO trough,, lama hatesta5demo 23melo trace tany
% SOs won't have the same length so we can get the poisions of troughs and
% then plot the raw signal [-2.5 +2.5] relative to the trough
idx = phenom.measures(:,1);
[~,temp] = cellfun(@(x) min(x),phenom.measures(:,2), 'Un',0);
trough_idx = cell2mat(cellfun(@(x,y) (x(y)) ,phenom.measures(:,1), temp, 'Un',0));


duration_to_cut = [-2.5 2.5];
pre_stim = duration_to_cut(1) * data.fsample;
post_stim = duration_to_cut(2) * data.fsample;

trough_idx(trough_idx<abs(pre_stim))=[]; % trimming the SO that came at the very beg. because we won't be able to see -2.5 before them and also from end
trough_idx(trough_idx+post_stim > length(data.trial{1}))=[];


% trlsx3 first col. beg. sample then end, offset: if zero then beg will be
% time0
trl = [pre_stim+trough_idx post_stim+trough_idx  repmat(pre_stim,length(trough_idx),1)];
cfg     = [];
cfg.trl = trl;
data_locked_trough = ft_redefinetrial(cfg,data);

% artf=ft_databrowser([],data_locked_trough); % visualize single trials

% erp of SO locked at trough
cfg             = [];
cfg.keeptrials  = 'yes';
data_locked_trough = ft_timelockanalysis(cfg, data_locked_trough);

mo_pretty_errorbar(data_locked_trough.time , squeeze(data_locked_trough.trial) , squeeze(data_locked_trough.trial).*0 ,  99);

% ylim([min(mean(data_locked_trough.trial,1)) max(mean(data_locked_trough.trial,1))])


plot(data_locked_trough.time, median(squeeze(data_locked_trough.trial),1))


%% manual locking
[~,temp] = cellfun(@(x) min(x),phenom.measures(:,2), 'Un',0);
trough_idx = cell2mat(cellfun(@(x,y) (x(y)) ,phenom.measures(:,1), temp, 'Un',0));

duration_samples = 1.5*data.fsample;
time_ax = -1.5:1/data.fsample:1.5;

trough_idx(trough_idx<duration_samples)=[];
trough_idx(trough_idx+duration_samples>length(data.trial{1}))=[];

event=[];
for i=1:length(trough_idx)
    trough = trough_idx(i);
    event(i,:) = data.trial{1}(trough-duration_samples : trough+duration_samples);
end



mo_pretty_errorbar(time_ax , event , (event).*0 ,  99);



%% SO spindle complexes (co-occurring of SO and spindles)
% SO
cfg=[]; result=[];
cfg.phenomenon = 'SO';
cfg.fun_seq = [{'band'} {'duration_zcrossings'} {'spec_thresholded'}];
cfg.data.fs = 200;
cfg.data = data;
[result.so_measures] = lv_detect_phenomenon(cfg);
%spindles
cfg=[];
cfg.phenomenon = 'spindle';
cfg.fun_seq = [{'band'}  {'spec_thresholded'}];
cfg.data.fs = 200;
cfg.data = data;
[result.spindle_measures] = lv_detect_phenomenon(cfg);

cfg=[];
cfg.phenomenon = 'SO_spindles complexes';
cfg.fun_seq = {'spec_thresholded'};
cfg.result=result;
result = lv_detect_phenomenon(cfg);

%% a testing for this and Active System Consolidation ASC testing is in PHASE_analysis_hilbert

%% actual coupling
% get SO and spindles (or any two phenom.) and then get the hilbert
% analytical signals to know the complex info during the time of the events
% then perform the actual coupling by taking the magnitude of the fast
% phenom. and the phase of the slower and put them in euler's form and then
% correct the result with surrogate data to have a modulation index which
% is complex valued: its magnitude reflects the coupling strength and its
% phase, the (slower phenom) phase for which there is the highest (faster phenom) amplitudes
% in this example we use SO as the slower phenom and spindles as the faster
% phenom.. if we want another phenom just consider the slower as so and the faster
% as spindle so you change the band of detecting the phenom and continue
% with the updated bands
% you may construct a vector of records were reactivation happened and
% check the coupling strength of that vs incorrect or control beaware of
% the channels used for SO and spindles ..
% SO
cfg=[];
cfg.phenomenon = 'SO';
cfg.fun_seq = [{'band'} {'duration_zcrossings'} {'spec_thresholded'}];
cfg.data.fs = 200;
cfg.data = data;
[so_measures] = lv_detect_phenomenon(cfg); % to get the timing of SO
%spindles
cfg=[];
cfg.phenomenon = 'spindle';
cfg.fun_seq = [{'band'}  {'spec_thresholded'}];
cfg.data.fs = 200;
cfg.data = data;
[spindle_measures] = lv_detect_phenomenon(cfg); % to get the timing of spindles

% SO_hilbert
cfg=[];
cfg.phenomenon = 'Instant analytical';
cfg.fun_seq = [{'band'}  {'spec_thresholded'}];
cfg.specs = [{'Instant analytical'} {'0.16 1.25'} {''} {'instant_phase_magnitude'} {'no'}]; % override the record with new values
cfg.data=data;
so_hilbert= lv_detect_phenomenon(cfg); % to get the phase of SO
% Spindle_hilbert
cfg=[];
cfg.phenomenon = 'Instant analytical';
cfg.fun_seq = [{'band'}  {'spec_thresholded'}];
cfg.specs = [{'Instant analytical'} {'12 16'} {''} {'instant_phase_magnitude'} {'no'}]; % override the record with new values
cfg.data=data;
spindle_hilbert = lv_detect_phenomenon(cfg); % to get the magnitude of spindles

% Coupling
cfg=[];
cfg.phenomenon = 'SO_spindles complexes';
cfg.fun_seq = {'coupled_SOspindle'};
cfg.data=data; % just to see data but with coupling we use the extracted measures
cfg.data.so_measures = so_measures;
cfg.data.spindle_measures = spindle_measures;
cfg.data.so_hilbert = so_hilbert;
cfg.data.spindle_hilbert= spindle_hilbert;
result = lv_detect_phenomenon(cfg); % results in: coupling_strength and phase_of_coupling



%% classification of motor imagery with CSP
% this is with one range mu range [8 12]HZ
warning('be careful of edges');
cfg_preprocessing                 = [];
cfg_preprocessing.bpfilter        = 'yes';
cfg_preprocessing.bpfreq          = [8 12];
data_bp= ft_preprocessing(cfg_preprocessing, data);

% csp build .. apply

% log_var manual .. lda


%% testing python classifier
load iris.dat

dat = iris(1:100,:); % this built-in data should be ready so no zscoring is done
dat((dat(:,end)==2),end)=-1;
X = dat(:,[1 3]);
Y = dat(:,end);
X = zscore(X, [], 1);
% visualising data
figure,
gscatter(X(:,1), X(:,2), Y ,'rgb','osd');
% classifier
sz = size(X);
id = randperm(sz(1));
TRAIN = X(id(1:round(sz(1)/2)),:);
TEST = X(id(round(sz(1)/2)+1:end),:);
TRAIN_GROUP = Y(id(1:round(sz(1)/2)));

TEST_GROUP = Y(id(round(sz(1)/2)+1:end));


classifier = py.importlib.import_module('py_classifiers');
prediction = classifier.NN(TEST,  TRAIN , TRAIN_GROUP, activation);

% lda
[outclass, err, posterior] = classify(TEST,  TRAIN , TRAIN_GROUP);
mean(outclass==TEST_GROUP)*100




%% Riemannian classification
% very nice based on Riemannian geometry and covariance matrices check the
% methods paper and the code in MotorImagery.m
% here we format our normal data into the needed format
% will reshape the covariance as a vector and put it in the place of
% channels and then reshape in the function for further analysis

% Riemannian geometry
cfg=[];
cfg.latency = [ 0 1.1]; im = ft_selectdata(cfg,im);
result_trn = im; result_tst=im;
result_trn.trial = im.trial(1:500,:,:); result_trn.trialinfo = clabel_hand_train(1:500);
result_tst.trial = im.trial(501:1000,:,:); result_tst.trialinfo = clabel_hand_train(501:1000); % tested .. just to test we split the trials this way

result_trn.trial = zscore(result_trn.trial,[],1);
result_tst.trial = zscore(result_tst.trial,[],1); % zscoring

% if we want to make the covariance on a time segment then it should be
% done here by adding other loop on time windows and the cov. is put instead of channels and the time points is as many windows as we use
% in the current case the whole time pts are used to get the cov. and so we
% will have a compressed time dim.  and the channel dim is the cov as a
% vector
temp=[];
for i=1:size(result_trn.trial,1), temp(i,:) = reshape(cov( squeeze(result_trn.trial(i,:,:))' ), 1,[]); end
result_trn.trial = temp; temp=[];
for i=1:size(result_tst.trial,1), temp(i,:) = reshape(cov( squeeze(result_tst.trial(i,:,:))' ), 1,[]); end
result_tst.trial = temp; temp=[];

cfg=[];
cfg.method = 'timextime';
cfg.classifier_type = {'Riemannian','tangent'}; % classifier, riemannian method
cfg.trn = result_trn; %has labels in trialinfo
cfg.tst = result_tst; cfg.folds=nan; cfg.do_parallel=0;
acc_mu{nn,1} = lv_classify(cfg);



%% toi classifiers, a classifier at every time point weighted by its CV performance and
% this temporal ensemble of classifiers is used to evaluate one timept from tst
% we train by firstly CV and getting the weights and then we train again
% to get the models then we apply on the testing set
% goes inside nn loop on sbjs
result_trn.trial = single(result_trn.trial); result_tst.trial = single(result_tst.trial);
% get the weights
cfg=[];
cfg.method = 'axtime';
cfg.classifier_type = {'lda'};
cfg.trn.trial = result_trn.trial; cfg.trn.trialinfo=result_trn.trialinfo(:,1);
cfg.folds=5; cfg.do_parallel=0;
classification_resultw(nn,:) = (squeeze(lv_classify(cfg)))';
toi_models.lv_weights(nn,:) = classification_resultw(nn,:) - 0.5; % deviation from chance is the weight
toi_models.lv_weights(nn,:) = toi_models.lv_weights(nn,:) ./ sum(toi_models.lv_weights(nn,:)); % softmax so that the weights sum to 1

% the bad classifiers are not just having the eopposite effect they could
% be totally random!! so you should pick something like 60% and only get
% the weights absed on these that you think carry the effect !! and use
% them with probabilities (0-1) or not depends !!
% % toi_models.lv_weights(nn,:) = classification_resultw(nn,:) - 0.6;
% % toi_models.lv_weights(nn,toi_models.lv_weights(nn,:)<0)=0;
% % toi_models.lv_weights(nn,:) = toi_models.lv_weights(nn,:) ./ sum(toi_models.lv_weights(nn,:)); % softmax so that the weights sum to 1

% get the models, the weights based on CV but the models are when we train with all trials
cfg=[];
cfg.method = 'axtime';
cfg.classifier_type = {'toi_classifiers','trn'};
cfg.trn.trial = result_trn.trial; cfg.trn.trialinfo=result_trn.trialinfo(:,1);
cfg.tst.trial = result_trn.trial.*0; cfg.tst.trialinfo=result_trn.trialinfo(:,1);
cfg.folds=nan; cfg.do_parallel=0;
temp = lv_classify(cfg); % 1_timeModels
toi_models.mdl(nn,:) = temp;
% get the results on testing set
cfg=[];
cfg.method = 'timextime'; % it's a curve but we use txt because we now got the models so we can evaluate once on the tst set like txt instead of looping on tst time
cfg.classifier_type =  {'toi_classifiers','tst',toi_models};
cfg.trn = result_trn; cfg.trn.trialinfo=result_trn.trialinfo(:,1);
cfg.tst = result_tst; cfg.tst.trialinfo=result_tst.trialinfo(:,1); cfg.folds=nan; cfg.do_parallel=0;
classification_result{nn,1} = lv_classify(cfg); %gives error be careful reduce data and don't use parallel ! and maybe use talldouble !


figure, plot( result_tst.time  ,classification_result{nn,1}); ylabel('auc'); xlabel('sleep time'); set(gca,'YDir','normal');
title(['classification for sbj: ' num2str(sbj(nn))])


%% positional embedding (FV_time)
% so it can be channel_time the positions of time features will encoded
% in the feature values not just ordered but actually encoded in the features values

% change the variables according to the dimensions of the data and then add the PE to the FVs you have
% i is the FV iterator 0:rows, pos is the time pt ... and d is the FV length
[ver, exec, loaded]	= pyversion;
%cd C:\Users\Mahmoud Eed\AppData\Local\Programs\Python\Python36\Scripts
% cd /d E:\Python\installed\Lib\site-packages
version = ('E:\Python\installed\python.exe');
pyversion(version) % use : pyversion(version) .. if you want to change to another version


% to produce pos_d positional embeddings
pos = 50; d = 128;
PE_fun = py.importlib.import_module('positional_encoding');
PE = PE_fun.positional_encoding(pos,d);

ar = double(py.array.array('d',PE.flatten('F').tolist()));
PE = reshape(ar,pos,d);

imagesc(PE)
xlabel('i'), ylabel('pos')
% then take this pos_FV and add the values to the feature values at
% different positions
% then we can take this d_pos and do a posxpos classification were the d
% are now the features ...



% Run after you edit py script: this will get the current module number and
% then rename it by incermenting the module number so that we can run the
% new one directly without restarting matlab
module_name = 'py_classifier_test';
my_library_file_path=['E:\Python\installed\Lib\site-packages\' module_name];
cd 'E:\Python\installed\Lib\site-packages';
file = dir(['*' module_name '*.py']); % to see the last file no.
current_file_no = str2num(string(extractBetween(file.name,module_name,'.py')));
my_library_file_path2 = [my_library_file_path num2str(current_file_no)];
movefile([my_library_file_path2 '.py'], ...
    [my_library_file_path num2str(current_file_no+1) '.py'] );
file_to_run = [module_name num2str(current_file_no+1)];



load fisheriris
activation_function = 'softmax';
idx = randperm(100); species = [ones(1,50) zeros(1,50)]; % got to be zeros and ones
TRAIN = meas(idx(1:50),:); TRAIN_GROUP = species(idx(1:50));
temp=TRAIN;
TEST = meas(idx(51:100),:); TEST_GROUP = species(idx(51:100));
TRAIN = py.numpy.array(TRAIN(:).');
TEST = py.numpy.array(TEST(:).');

d_fun = py.importlib.import_module(file_to_run);
prediction = d_fun.NN(TEST, TRAIN, TRAIN_GROUP, activation_function);

estimated = double(py.array.array('d',prediction.flatten('F').tolist()));
mean(estimated==TEST_GROUP)*100



% matlab to python
npA = py.numpy.array(A(:).');
% python to matlab
ar = double(py.array.array('d',npA.flatten('F').tolist()));
sss = reshape(ar,3,3);

% quick test
d_fun = py.importlib.import_module('test_tf');
xx = d_fun.dd;

% lda
[outclass, err, posterior] = classify(TEST,  TRAIN , TRAIN_GROUP);
mean(outclass==TEST_GROUP')*100

% to test in python commands because matlab needs  to restart everytime we edit py code
% import py_classifier_test as fun
% fun.NN(0,[[1,2,3],[4,5,6]], 0,0)



%% signals alignment
% calculate crosscorrelation
a = x; b = x;
c = [ zeros(1,length(a)-1) a ];
d = [ b zeros(1,length(a)-length(b)+length(a)-1) ];
e = fft(c);
f = fft(d);
g = e.*conj(f);
h = ifft(g);
% or
xcorr(x,x)

% shifted signals example
% from example here: https://www.mathworks.com/help/signal/ug/compensate-for-delay-and-distortion-introduced-by-filters.html
% phase_shifted_signals.trial = [xn;xf];
% phase_shifted_signals.time = tn;
% save phase_shifted_signals phase_shifted_signals

load phase_shifted_signals phase_shifted_signals
% phase_shifted_signals.trial(1,:) = -phase_shifted_signals.trial(1,:); %
% to see the direction effect .. dife
dat.trial(:,1,:) = phase_shifted_signals.trial;
dat.trial(:,2,:) = phase_shifted_signals.trial;
dat.trial(:,3,:) = phase_shifted_signals.trial;
% dat.trial = [dat.trial ; dat.trial];
dat.time = phase_shifted_signals.time;
% with lv
cfg=[];
cfg.data = dat; % data.time all time extending beyond trial limit
cfg.roi =[0.2 0.6]; cfg.shift=170;  cfg.method='topos';
dat2 = lv_align(cfg);

figure,
plot(dat2.time,squeeze(dat2.trial(1,3,:)))
hold on, plot(dat2.time,squeeze(dat2.trial(2,3,:)))
figure,
plot(dat.time,squeeze(dat.trial(1,3,:)))
hold on, plot(dat.time,squeeze(dat.trial(2,3,:)))

% CORAL test (alignment of features in the feature space not the time
% points)
% didn't see the difference because these signals are phase shifted and when CORAL scatter the features
% it won't care about the temporal information so eventhough they are phase
% shifted they will still occupy the same feature space so it won't change
% the signal so we need to change values not phase .. now we can think of
% EA and CORAL as magnitude alignment and xcorrelation as phase alignment
cfg=[];  rng(1)
cfg.trn.trial=repmat(phase_shifted_signals.trial(1,:),1000,1) + randn(1000,size(phase_shifted_signals.trial,2))./50;
cfg.tst.trial=repmat(phase_shifted_signals.trial(2,:)/2,1000,1) + randn(1000,size(phase_shifted_signals.trial,2))./50;
cfg.method='CORAL';
features_data = lv_feature_extractor(cfg)
plot(dat.time,cfg.trn.trial)
hold on, plot(phase_shifted_signals.time,cfg.tst.trial)

figure,
plot(phase_shifted_signals.time,squeeze(features_data.trial))
hold on, plot(phase_shifted_signals.time,cfg.tst.trial)


%% align new data with lv xcorrelation
% general params
start_up
addDirs
fstruct = dir([rawdir '/' 'part*.eeg']); h = struct2cell(fstruct);
sbj = unique(str2double( (cellfun(@(x) (strtok(x,'part_')), (h(1,:))','Un', 0))' ));
sbj(ismember(sbj,[ 2 19 13 25,3     4     5     8    17    20    23    24    26]))=[]; % bad sbj to remove, sbj2 had different reference and channels..
cleaned_path = [preprocdir 'final_cleaned_after_inspection\part'];
almostraw = 'D:\sul''s code\Matt\sleep\erps\Organised\New exp/RawData/allSleep/data/cleaned/';

for nn=1:numel(sbj)
    fprintf(['\n lv: Working on pipeline for sbj: ' num2str(sbj(nn)) ' \n']);
    im = lv_load([cleaned_path num2str(sbj(nn)) '_img_manual_cleaned_N0'],'trial');
    sleep = lv_load([cleaned_path num2str(sbj(nn)) '_sleep_manual_cleaned_N3'],'trial');
    shift = 50; % 50ms maxshift both directions 50ms right and 50ms left
    % data trimming
    cfg=[]; cfg.latency=[0-(shift/1000) 1.1+(shift/1000)]; % shift is considered
    im=ft_selectdata(cfg,im);
    cfg.latency=[0-(shift/1000) 2.5+(shift/1000)];
    sleep=ft_selectdata(cfg,sleep);
    % aligning
    cfg=[]; cfg.roi=[0 1.1]; cfg.data=im;
    cfg.shift=shift;
    im = lv_align(cfg);
    cfg.roi=[0 2.5]; cfg.data=sleep;
    sleep = lv_align(cfg);
    
    lv_save([cleaned_path num2str(sbj(nn)) '_img_manual_cleaned_aligned_N0'],im,'trial');
    lv_save([cleaned_path num2str(sbj(nn)) '_sleep_manual_cleaned_aligned_N3'],sleep,'trial');
end


%% testing RNN in matlab
load HumanActivityTrain
XTrain


X = XTrain{1}(1,:);
classes = categories(YTrain{1});

figure
for j = 1:numel(classes)
    label = classes(j);
    idx = find(YTrain{1} == label);
    hold on
    plot(idx,X(idx))
end
hold off

xlabel("Time Step")
ylabel("Acceleration")
title("Training Sequence 1, Feature 1")
legend(classes,'Location','northwest')


% dividing into trials ...
% X = XTrain{1}(:,1:64400); % from one sbj as training and many cells as trials ...
% Y = YTrain{1}(:,1:64400);
% XX = mat2cell(X,3,repmat(10,1,6440));
% YY = mat2cell(Y,1,repmat(10,1,6440));


numFeatures = 3;
numHiddenUnits = 200;
numClasses = 5;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits,'OutputMode','sequence')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];


options = trainingOptions('adam', ...
    'MaxEpochs',60, ...
    'GradientThreshold',2, ...
    'Verbose',0, ...
    'Plots','training-progress');


net = trainNetwork(XX,YY,layers,options);  % XTrain: guess here it was sbjs for us it will be trials

%test
load HumanActivityTest
figure
plot(XTest')
xlabel("Time Step")
legend("Feature " + (1:numFeatures))
title("Test Data")
YPred = classify(net,XTest);
acc = sum(YPred == YTest)./numel(YTest)

figure
plot(YPred,'.-')
hold on
plot(YTest)
hold off

xlabel("Time Step")
ylabel("Activity")
title("Predicted Activities")
legend(["Predicted" "Test Data"])



% on EEG
for i=1:size(data.trial,3), temp{1,i} = squeeze(data.trial(i,:,:));
    labels{1,i} = repmat(GROUP_TRAIN(i),1,size(data.trial(i,:,:),3)); end
data.trial = temp; clear temp;
% RNN train
numFeatures = size(data.trial,2);
numHiddenUnits = 200;
numClasses = length(unique(GROUP_TRAIN));
layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits,'OutputMode','sequence')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];
options = trainingOptions('adam', ...
    'MaxEpochs',10, ... % was 60
    'GradientThreshold',2, ...
    'Verbose',0, ...
    'Plots','training-progress');
net = trainNetwork(data.trial,labels,layers,options);  % XTrain: guess here it was sbjs for us it will be trials
% RNN test
YPred = classify(net,TEST); % test shouyld be in the same format
% YPred  is now across time ...

%% lv_align with a reference which makes sense because many trials will result in near zero shift which is noisy
% lv_align choosing a reference signal (classless it should be a global reference
cfg=[]; cfg.latency = [ 0 1.1]; data = ft_selectdata(cfg,imHold); coeff=[];
progressbar = ParforProgressbar(size(data.trial,1),'title', 'trials'' correlations progress:');
parfor i=1:size(data.trial,1)
    for j=i:size(data.trial,1)
        c = normxcorr2( squeeze(data.trial(i,:,:)), squeeze(data.trial(j,:,:)) );
        
        % only looking at the logical parts without edge artifacts of the 2D convolution, so we remove the padded convolution results
        c(:,[1:size(data.trial,3)-1 end-(size(data.trial,3)-1)+1:end])=[];
        c([1:size(data.trial,2)-1 end-(size(data.trial,2)-1)+1:end],:)=[];
        coeff = [coeff   abs(c)]; % see the commented part above we just do it in a row for parfor to work
        %        lag(i,j) = -lag(j,i);
    end
    progressbar.increment();
end
delete(progressbar);
count=0; % refilling with the results from parfor
for i=1:size(data.trial,1)
    for j=i:size(data.trial,1)
        count=count+1; coeff_all(j,i) = coeff(count);
        coeff_all(i,j) = coeff_all(j,i);
    end
end

coeff = mean(coeff_all,2);

[val,ref_trl(nn,1)] = max(coeff)
% aligning with lv with xcorrelation
cfg=[];
cfg.data = imHold; % data.time all time extending beyond trial limit
cfg.roi =[0 1.1]; cfg.method='topos';
cfg.shift = 200;
cfg.ref_trl = squeeze(data.trial(ref_trl(nn,1),:,:));
im_aligned = lv_align(cfg);

figure,
subplot(121)
imagesc(im_aligned.time,1:size(im_aligned,1),[squeeze(im_aligned.trial(clabel_hand_train==1,4,:)) ; squeeze(im_aligned.trial(clabel_hand_train==2,4,:))]); % pick a channel and see trlxtime
caxis([-20 20])
subplot(122),
cfg=[]; cfg.latency = [ 0 1.1]; imHold2 = ft_selectdata(cfg,imHold);
imagesc(imHold2.time,1:size(imHold2,1),[squeeze(imHold2.trial(clabel_hand_train==1,4,:)) ; squeeze(imHold2.trial(clabel_hand_train==2,4,:))]);
caxis([-20 20])

%% hilbert extracts based on cosine phase
fs=100; t = 0:1/fs:1;
signal = 4*sin(2*pi*5*t);
plot(t,signal)
hold on,
plot(t,angle(hilbert(signal')));
figure,%power
plot(t,abs(hilbert(signal')).^2);



%% ITPC on small dataset and compare to mike's lecture 150
% formatting data to do the tf analysis
load('sampleEEGdata.mat'); % good to use this to compare to his analyses
data = [];
data.trial = permute(EEG.data,[3 1 2]);
data.trialinfo=ones(size(data.trial,1),1);
temp = 1:length(EEG.times):length(EEG.times)*size(data.trial,1);
data.sampleinfo = [temp' temp'+length(EEG.times)-1];
data.time = EEG.times/1000;
data.method = 'itpc';
data.label = arrayfun(@(x) (x.labels), EEG.chanlocs, 'Un',0); data.dimord='rpt_chan_time';

load('sampleEEGdata.mat'); % good to use this to compare to his analyses
data = [];
data.trial = permute(EEG.data,[3 1 2]);
data.time = EEG.times/1000;
data.method = 'itpc';
data.label = arrayfun(@(x) (x.labels), EEG.chanlocs, 'Un',0); data.dimord='rpt_chan_time';
cfg=[]; cfg.data=data;
data = lv_to_ftstruct(cfg);


TF_struct = lv_tf(data, 0, 0); %data, do_stats, do_plot
% group lvl TF   .. use on one channel for one figure
frequencies = [0.5 25];
TF_struct.freq = linspace(frequencies(1),frequencies(end),2*(1+frequencies(end)-frequencies(1))); % same as inside lv_tf
TF_struct.label = data.label;
TF_struct.time = data.time;
% lv_tf(TF_struct, 0, 1); %data, do_stats, do_plot ... it's working but you need mike's position of channels and it's fine not to do it because he works on one channel

% on one channel like mike
channel2use = 27; %o1
itpc = squeeze(TF_struct.powspctrm); figure,
contourf(TF_struct.time,TF_struct.freq,squeeze(itpc(channel2use,:,:)),50,'linecolor','none') % to smooth the result
set(gca,'clim',[0 .7],'xlim',[-.5 1.3])
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
title([ 'ITPC' ])

% taking slope
fs=100; t = 0:1/fs:1;
signal = 4*cos(2*pi*2*t);
plot(t,signal)
figure,
plot(t(1):1/fs:t(end-1),diff(signal)) %y2-y1 but gradient is good to get the same signal's length
figure,
plot(t,gradient(signal))

%% lateralised data and high,low variance data and GED
fs=1000; t=0:1/fs:1;
s = 5*sin(2*pi*50*t);
plot(t,s)
channels = 6;

% high variance signal
r = rand(channels,length(s)); % in the range [0 1]
r(:,700:800)=r(:,700:800)+2;
r(1:3,710:750)=r(1:3,710:750)+0.25;
s1 = bsxfun(@times,r,s);

% hold on, plot(t,-abs(hilbert(s1(1,:)'))')
% low variance signal
r = rand(channels,length(s)); % in the range [0 1]
r(:,700:800)=r(:,700:800)-0.5;
r(1:3,760:860)=r(1:3,760:860)-0.25;
s2 = bsxfun(@times,r,s);

% changing the values of some channels to simulate data and increase the rank
s1(5,:) = circshift(s1(5,:),2); s1(6,:) = circshift(s1(6,:),10);
s2(5,:) = circshift(s2(5,:),8); s2(6,:) = circshift(s2(6,:),15);

% re-run for both
covs1 = cov(s1'); %timexchannel
imagesc(covs1) % high variance with high values gives high covariance

figure
covs2 = cov(s2'); %timexchannel
imagesc(covs2)% low variance with low values gives low covariance


[W,val] = eig(covs1,covs2); % generalized decomposition
imagesc(real(W))
figure
imagesc(real(val))
[~,id] = sort(diag(val),'descend');


% GED separated dims
s1_projected = [W(:,id(1))' * s1 ; W(:,id(2))' * s1]; % using two eigvectors that maximally separate the data
s2_projected = [W(:,id(1))' * s2 ; W(:,id(2))' * s2];
figure
scatter(t,s1_projected(1,:),'b')
hold on, scatter(t,s2_projected(1,:),'r')
figure
scatter(t,s1(1,:),'b')
hold on, scatter(t,s2(1,:),'r')

s1_projected=s1_projected.^2;
s2_projected=s2_projected.^2;

% lateralised:
% low variance signal
r = rand(channels,length(s)); % in the range [0 1]
r(:,700:800)=r(:,700:800)+3.5;
r(1:3,760:860)=r(1:3,760:860)+0.25;
s2 = bsxfun(@times,r,s);
% changing the values of some channels to simulate data and increase the rank
s2(5,:) = circshift(s2(5,:),8); s2(6,:) = circshift(s2(6,:),15);

s1 = abs(hilbert(s1'))'; % lateralised effect (positive from one channel
s2 = -abs(hilbert(s2'))'; % lateralised effect (positive from one channel
figure,
plot(t,s2)
hold on, plot(t,s1)
% now run the re-run part to do the same but foe lateralised signals


% projected features (2d) two channels
figure
scatter(s1_projected(1,:),s1_projected(2,:),'b')
hold on, scatter(s2_projected(1,:),s2_projected(2,:),'r')
scatter(s1(1,:),s1(2,:),'b') % original
hold on, scatter(s2(1,:),s2(2,:),'r')

% amplitude to variance .. function
sign_vec = mod(1:size(s1,2), 2); sign_vec(sign_vec==0)=-1;
s1(s1<0) = s1(s1<0).^-1;
s2(s2<0) = s2(s2<0).^-1;

s1 = bsxfun(@times,s1,sign_vec);
s2 = bsxfun(@times,s2,sign_vec);

% inv of the function
s1_projected = bsxfun(@times,s1_projected,sign_vec);
s2_projected = bsxfun(@times,s2_projected,sign_vec);

s1_projected(s1_projected<0) = s1_projected(s1_projected<0).^-1;
s2_projected(s2_projected<0) = s2_projected(s2_projected<0).^-1;

% load fisheriris
% Mdl = fitcdiscr(meas,species)

% doing GED for between, within because this should have different cov

all_covs = cov([s1' s2'],1 ); %timexchannel the ,1 adds normalisation by the no. timpts that will be multiplied and removed later
% from that big cov we take the within and between class
figure, imagesc(all_covs)

Swithin = ((size(s1,2).*all_covs(1:channels,1:channels)) + (size(s2,2).*all_covs(channels+1:end,channels+1:end))); % for same number of timepts for different number we will multiply by n-1 where n is the no. timepts
figure, imagesc(Swithin)
% Sbetween=[];
% for i=1:channels
%     for j=i:channels
%         temp = cov([s1(i,:)' s2(j,:)']);
%         Sbetween(i,j) = temp(1,2); %all_covs(i,j+channels);
%         Sbetween(j,i) = Sbetween(i,j);
%     end
% end
% Sbetween=size(s1,2).*Sbetween;
Sbetween = cov([s1';s2']);
Sbetween=size(s1,2).*Sbetween;
figure, imagesc(Sbetween)

[W,val] = eig(Sbetween,Swithin); % generalized decomposition
[~,id] = sort(diag(val),'descend');

% use % GED separated dims ... to project and plot


%% with matthias' ldaBeamformer
% To obtain the LDA beamformer, we need the spatial pattern P and the
% corresponding continuous or epoched data (NOT the averaged data).
erp_pattern = mean(s1-s2, 2); % difference of classes erps is considered instead of cov_between .. and compressed in time
epoched_data(:,:,1) = s1; epoched_data(:,:,2) = s2;% 3D matrix (channels x time samples x epochs)
[w1,lda1,C1] = LDAbeamformer(erp_pattern,epoched_data);


figure
scatter(t,s1(1,:),'b')
hold on, scatter(t,s2(1,:),'r')
figure
scatter(t,lda1(1,:),'b')
hold on, scatter(t,lda1(2,:),'r')

%% projecting on actual data
% GED projection to seperate classes
channels=size(im.trial,2);
t=im.time(nearest(im.time,0.5):nearest(im.time,1.1));
s1 = squeeze(mean(im.trial(clabel_hand_train==1,:,nearest(t,0.5):nearest(t,1.1)),1));
s2 = squeeze(mean(im.trial(clabel_hand_train==2,:,nearest(t,0.5):nearest(t,1.1)),1));
all_covs = cov([s1' s2'],1 ); %timexchannel the ,1 adds normalisation by the no. timpts that will be multiplied and removed later
% from that big cov we take the within and between class
%     figure, imagesc(all_covs)

Swithin = ((size(s1,2).*all_covs(1:channels,1:channels)) + (size(s2,2).*all_covs(channels+1:end,channels+1:end))); % for same number of timepts for different number we will multiply by n-1 where n is the no. timepts
%     figure, imagesc(Swithin)

% Sbetween=[]; laaa 3ayzeno classless !! la2n 3adad el observtion elel hwa hena el time mmkn yb2a mo5talef aslan fa el mafrod classless
% for i=1:channels
%     for j=i:channels
%         temp = cov([s1(i,:)' s2(j,:)']);
%         Sbetween(i,j) = temp(1,2); %all_covs(i,j+channels);
%         Sbetween(j,i) = Sbetween(i,j);
%     end
% end
% Sbetween=size(s1,2).*Sbetween;
Sbetween = cov([s1';s2']);
Sbetween=size(s1,2).*Sbetween;

%     figure, imagesc(Sbetween)

[W,val] = eig(Sbetween,Swithin); % generalized decomposition
[~,id] = sort(diag(val),'descend');

new_im=[];
for i=1:size(im.trial,1)
    new_im.trial(i,:,:) = W(:,id(1))' * squeeze(im.trial(i,:,:));
end

cfg = [];
cfg.classifier      = 'lda';
cfg.metric          = 'accuracy';   cfg.preprocess      = {'zscore'};
cfg.hyperparameter           = [];
cfg.hyperparameter.lambda    = 'auto';
rng(1)
[txtAuc(nn,:,:),result{nn,1} ] = mv_classify_timextime(cfg,im.trial  , clabel_hand_train' );
[txtAuc_GED(nn,:,:),result{nn,1} ] = mv_classify_timextime(cfg,new_im.trial  , clabel_hand_train' );

%% domain alignment with subspace alignment based on PCs align
% 2d data
mu = [0 0];
Sigma = [1 1.5; 1.5 3]; rng(1)
trn_c1 = mvnrnd(mu,Sigma,100); %outlier = find(sum(trn_c1,2)==max(sum(trn_c1,2),[],1)); trn_c1(outlier,:)=trn_c1(outlier,:)+[10 90];
trn_c2 = trn_c1-5;

trn = [trn_c1;trn_c2];
trn = zscore(trn,[],1);
trn = trn - repmat(mean(trn,1),size(trn,1),1);
trn_c1 = trn(1:100,:); trn_c2 = trn(101:end,:);

figure, plot(trn_c1(:,1),trn_c1(:,2),'+')
hold on, plot(trn_c2(:,1),trn_c2(:,2),'o')

rotation = -60; % change rotation to test new combinations .. should be negative for  the features to take the same side for ddifferent classes 
%o.w. the effect means that the classes are swapped and there is no classless rotation that will solve that and this is a problem with the data itself
T = [cosd(rotation) -sind(rotation);sind(rotation) cosd(rotation)];
tst = T*trn';
tst_c1 = tst(:,1:100)'; tst_c2 = tst(:,101:end)';

tst = [tst_c1;tst_c2];
tst = zscore(tst,[],1);
tst_c1 = tst(1:100,:); tst_c2 = tst(101:end,:);

plot(tst_c1(:,1),tst_c1(:,2),'+')
hold on, plot(tst_c2(:,1),tst_c2(:,2),'o')

% pca
[eigVects_target] = pca(trn); [eigVects_source] = pca(tst);
trn = (eigVects_target(:,1)' * trn')';
tst = (eigVects_target(:,1)' * tst')';
figure, plot(trn(1:100,1),zeros(1,100),'+')
hold on, plot(trn(101:end,1),zeros(1,100),'o')
% CCA
[a,b] = canoncorr(trn,tst); % find the linear trnasformation a that transforms features of trn
% and b that transforms features of tst to maximise the correlation between
% pairs of corresponding features 
trn = trn*a;
tst = tst*b;
figure, plot(trn(1:100,1),trn(1:100,2),'+')
hold on, plot(trn(101:end,1),trn(101:end,2),'o')
hold on, plot(tst(1:100,1),tst(1:100,2),'+')
plot(tst(101:end,1),tst(101:end,2),'o')

% SA align
cfg=[];
cfg.trn.trial=trn; cfg.tst.trial=tst;
cfg.method = 'CORAL';
transformed = lv_feature_extractor(cfg);

trn_c1_t = transformed.source(1:100,:); trn_c2_t = transformed.source(101:end,:);
tst_c1_t = transformed.target(1:100,:); tst_c2_t = transformed.target(101:end,:);

figure, plot(trn_c1_t(:,1),trn_c1_t(:,2),'+')
hold on, plot(trn_c2_t(:,1),trn_c2_t(:,2),'o')
plot(tst_c1_t(:,1),tst_c1_t(:,2),'+')
hold on, plot(tst_c2_t(:,1),tst_c2_t(:,2),'o')

% from domain_align toolbox .. found that data isn't centered inside the
% function
ftAll = [trn;tst]; maSrc = 0;%logical([ones(size(trn,1),1);zeros(size(tst,1),1)]);
maLabeled = 0;%ones(size(ftAll,1),1);
target = 0;%[ones(100,1);2*ones(100,1);ones(100,1);2*ones(100,1)];

[ftAllNew,transMdl] = ftTrans_ssa(ftAll,maSrc,target,maLabeled); % you can try others like ftTrans_tca

trn_c1_t = ftAllNew(1:100,:); trn_c2_t = ftAllNew(101:200,:);
tst_c1_t = ftAllNew(201:300,:); tst_c2_t = ftAllNew(301:end,:);

figure, plot(trn_c1_t(:,1) ,'+')
hold on, plot(trn_c2_t(:,1) ,'o')
plot(tst_c1_t(:,1) ,'+')
hold on, plot(tst_c2_t(:,1) ,'o')

% CORAL
cfg=[];
cfg.trn.trial=trn; cfg.tst.trial=tst;
cfg.method = 'CORAL_manual'; % CORAL
transformed = lv_feature_extractor(cfg);

trn_c1_t = transformed.trial(1:100,:); trn_c2_t = transformed.trial(101:end,:);
figure, plot(trn_c1_t(:,1),trn_c1_t(:,2),'+')
hold on, plot(trn_c2_t(:,1),trn_c2_t(:,2),'o')
plot(tst_c1(:,1),tst_c1(:,2),'+')
hold on, plot(tst_c2(:,1),tst_c2(:,2),'o')

% % CSP testing on the scattered points just to vissualise its effect and it was good
% % because we want to see the effect on scattered points we will consider
% % this as one trial from each class and the old trials points as the time points in the trial
% cfg=[]; % building CSP filters
% cfg.method='csp';
% cfg.step='calculate';
% cfg.data.trial(1,:,:)=trn'; cfg.data.trial(2,:,:)=tst';
% cfg.data.trialinfo=[1;2];
% comp = lv_component_analysis(cfg);
% comp.id = comp.id(1:2); % choosing csp filters from sorted
% cfg.step='transform'; % applying the chosen filters to transform data,, in a CV paradigm when we use all data
% % to build the csp filters there will be leakage but we don't care much about it here
% cfg.comp=comp;
% ftAllNew = lv_component_analysis(cfg);
% trn_c1_t = squeeze(ftAllNew(1,:,1:100))'; trn_c2_t = squeeze(ftAllNew(1,:,101:end))';
% tst_c1_t = squeeze(ftAllNew(2,:,1:100))'; tst_c2_t = squeeze(ftAllNew(2,:,101:end))';
% figure, plot(trn_c1_t(:,1),trn_c1_t(:,2) ,'+')
% hold on, plot(trn_c2_t(:,1),trn_c2_t(:,2) ,'o')
% plot(tst_c1_t(:,1),tst_c1_t(:,2) ,'+')
% hold on, plot(tst_c2_t(:,1),tst_c2_t(:,2) ,'o')

%% von-mises distribution for handling circular data
% generate random angles
kappa=5; theta=0; % kappa is the concentration and is like the inverse of the varianceand theta is the preferred angle (mean)
alpha1 = circ_vmrnd(theta, kappa, 100);
kappa = 5; theta=pi;
alpha2 = circ_vmrnd(theta, kappa, 100);

figure, circ_plot(alpha1,'pretty','ro',true,'linewidth',2,'color','r') % class1 is red
hold on, circ_plot(alpha2,'pretty','bo',true,'linewidth',2,'color','b')
title('original classes')

% alphax = circ_vmrnd(theta, 0, 100);
% circ_plot(alphax,'pretty','bo',true,'linewidth',2,'color','b')
% [pval z] = circ_rtest(alphax)

[mu1] = circ_mean(alpha1); [mu2] = circ_mean(alpha2);
kappa1 = circ_kappa(alpha1); kappa2 = circ_kappa(alpha2);

% new samples to classify
kappa = 0.1; theta=pi;
alpha_tst = circ_vmrnd(theta, kappa, 400);
figure, circ_plot(alpha_tst,'pretty','ro',true,'linewidth',2,'color','r')
title('testing samples')

likelihood1 = circ_vmpdf(alpha_tst, mu1, kappa1); % directly from the toolbox but it's univariate
likelihood2 = circ_vmpdf(alpha_tst, mu2, kappa2);

posterior = [likelihood1./(likelihood1+likelihood2) likelihood2./(likelihood1+likelihood2)];

est = 1*(posterior(:,1)>posterior(:,2)) + 2*(posterior(:,1)<posterior(:,2));
figure, circ_plot(alpha_tst(est==1),'pretty','ro',true,'linewidth',2,'color','r')
hold on, circ_plot(alpha_tst(est==2),'pretty','bo',true,'linewidth',2,'color','b')
title('classification output')
% evaluate pdf for multivariate.. trials .. it doesn't work so if we have
% many channels we can take the circular mean and then work with one feature
% % likelihood1 = exp([kappa1 kappa1]*cos([alpha_tst alpha_tst]-[mu1 mu1])')/(2*pi*besseli(0,mean([kappa1 kappa1])));
% trial multivariate (sum( exp(kappa1.*cos(alpha_tst)) )./(2*pi))*2*pi;
% fun = @(x) exp(kappa1*cos(alpha_tst));
% q = (integral(fun,0,2*pi,'ArrayValued',1)./(2*pi))*2*pi;
% % likelihood2 = exp([kappa2 kappa2]*cos([alpha_tst alpha_tst]-[mu2 mu2])')/(2*pi*besseli(0,mean([kappa2 kappa2])));
% % likelihood1=likelihood1'; likelihood2=likelihood2';
% with 2d lda 
[mu1] = mean([cos(alpha1) sin(alpha1)],1); [mu2] = mean([cos(alpha2) sin(alpha2)],1);
cov1 = cov([cos(alpha1) sin(alpha1)]); cov2 = cov([cos(alpha2) sin(alpha2)]);
likelihood1 = mvnpdf( [cos(alpha_tst) sin(alpha_tst)], mu1, cov1); % directly from the toolbox but it's univariate
likelihood2 = mvnpdf([cos(alpha_tst) sin(alpha_tst)], mu2, cov2);
posterior = [likelihood1./(likelihood1+likelihood2) likelihood2./(likelihood1+likelihood2)];
est = 1*(posterior(:,1)>posterior(:,2)) + 2*(posterior(:,1)<posterior(:,2));
figure, circ_plot(alpha_tst(est==1),'pretty','ro',true,'linewidth',2,'color','r')
hold on, circ_plot(alpha_tst(est==2),'pretty','bo',true,'linewidth',2,'color','b')
title('2d LDA classification output')

%% t-distribution lda not normal distribution 
% t-distribution is less sensitive to outliers so it may give good results 
% generate data
mu = [0 0];
Sigma = [1 1.5; 1.5 3]; rng(1)
trn_c1 = mvnrnd(mu,Sigma,100); %outlier = find(sum(trn_c1,2)==max(sum(trn_c1,2),[],1)); trn_c1(outlier,:)=trn_c1(outlier,:)+[10 90];
trn_c2 = trn_c1-5;
trn = [trn_c1;trn_c2];
trn = zscore(trn,[],1);
trn_c1 = trn(1:100,:); trn_c2 = trn(101:end,:);

figure, plot(trn_c1(:,1),trn_c1(:,2),'+')
hold on, plot(trn_c2(:,1),trn_c2(:,2),'o')

mu1 = mean(trn_c1,1); mu2 = mean(trn_c2,1); 
% new samples to classify
mu = [0 0];
Sigma = [4 1.5; 1.5 4]; rng(1)
tst = mvnrnd(mu,Sigma,1000); 
hold on, plot(tst(:,1),tst(:,2),'+')
title('testing samples')

df=10000; % degree of freedom .. equals gaussian when df=inf.. and cauchy distribution when df=1 .. a parameter to optimise
likelihood1 = mvtpdf(bsxfun(@minus,tst,mu1),cov(trn_c1),df); % we subtract the mean of the class because mvtpdf has the assumption that the t-dist. is at zero mean 
% and in equation 2.159 the mean was subtracted 
likelihood2 = mvtpdf(bsxfun(@minus,tst,mu2),cov(trn_c2),df);
posterior = [likelihood1./(likelihood1+likelihood2) likelihood2./(likelihood1+likelihood2)];
est = 1*(posterior(:,1)>posterior(:,2)) + 2*(posterior(:,1)<posterior(:,2));
figure, plot(tst(est==1,1),tst(est==1,2),'+','color','b')
hold on, plot(tst(est==2,1),tst(est==2,2),'o','color','r')
title('2d LDA classification output')
hold on, plot(trn_c1(:,1),trn_c1(:,2),'+') % original
hold on, plot(trn_c2(:,1),trn_c2(:,2),'o')
% testing outliers .. Fig 2.16
rng(1)
x = mvnrnd(0,1,300)'; x = [x repmat(8,1,10) repmat(7,1,10) repmat(11,1,10)];% adding outliers
figure,subplot(121),
h=histogram(x,'Normalization','pdf'); sum(h.Values) % the prob. not the main thing we want to see the shift
xlim([-5 10])
% normal fit
[m,s] = normfit(x);
y = normpdf(x,m,s);
hold on, plot(x,y,'.'); title('normal dist. fit')
% t-fit
subplot(122),
h=histogram(x,'Normalization','pdf'); sum(h.Values) % the prob. not the main thing we want to see the shift
xlim([-5 10])
y = tpdf(x-mean(x),100000000000);
hold on, plot(x,y,'.');  title('t dist. fit') 

%% whiten data with PCA and there are different methods like cholesky method 
mu = [0 0];
Sigma = [1 1.5; 1.5 3]; rng(1)
trn_c1 = mvnrnd(mu,Sigma,100); %outlier = find(sum(trn_c1,2)==max(sum(trn_c1,2),[],1)); trn_c1(outlier,:)=trn_c1(outlier,:)+[10 90];
 
trn = zscore(trn_c1,[],1); 

figure, subplot(3,1,1); plot(trn(:,1),trn(:,2),'+'), title('original');
xlim([-3 3]), ylim([-3 3])

[ eigVects,~, eigVals,~] = pca(trn);

trn2 = (eigVects' * trn')'; % decorrelate data by transforming with the PCs
subplot(3,1,2); plot(trn2(:,1),trn2(:,2),'+'), title('original');
trn2(:,1)=trn2(:,1)./sqrt(eigVals(1)); % dividing by the std to make it one such that 
% it becomes a circle so because the eigVal is the variance we take the
% sqrt and divide by it for the respective component..
trn2(:,2)=trn2(:,2)./sqrt(eigVals(2));
subplot(3,1,3); plot(trn2(:,1),trn2(:,2),'r+'), title('whitened');
xlim([-3 3]), ylim([-3 3])
% verify
var(trn2,[],1)
mean(trn2,1)
%% lv align
% lv_align choosing a reference signal (classless it should be a global reference
% applied to wake erps after beamformer projection using sleep spatial filter
coeff=[]; data = result_trn; % cut to region of interest(trial 1.5sec.)
% temp = permute(erp,[3 2 1]); %ppnt_time
% id=mod(1:size(erp,3), 2); data.trial=[];
% data.trial(:,1,:) = [squeeze(temp(:,:,1)) ; squeeze(temp(:,:,2))]; 
data.trial=[]; data.trial(:,1,:) = axtimeAcc; % for axtime classification

hold_full_signals = data;
data.trial = data.trial(:,:,nearest(data.time,0):nearest(data.time,1.1)); % region of interest

for i=1:size(data.trial,1)
    for j=i:size(data.trial,1)
        c = reshape( normxcorr2( squeeze(data.trial(i,:,:)), squeeze(data.trial(j,:,:)) ),size(data.trial,2),[] );
        % only looking at the logical parts without edge artifacts of the 2D convolution, so we remove the padded convolution results
        c(:,[1:size(data.trial,3)-1 end-(size(data.trial,3)-1)+1:end])=[];
        c([1:size(data.trial,2)-1 end-(size(data.trial,2)-1)+1:end],:)=[];
        coeff = [coeff   abs(c)]; % see the commented part above we just do it in a row for parfor to work
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
cfg=[];
cfg.data = hold_full_signals; % data.time all time extending beyond trial limit
cfg.roi =[0 1.1]; cfg.method='topos';
cfg.shift = 400; cfg.classful=1;
cfg.ref_trl = reshape(squeeze(data.trial(ref_trl,:,:)),size(data.trial,2),[]);
erp = lv_align(cfg);
lv_pretty_errorbar(erp.time, squeeze(erp.trial(1:length(sbj),:,:)), ...
    squeeze(erp.trial(length(sbj)+1:end,:,:)), 1);

lv_vec2corr(w1, w2, 'xaxxis','yaxis' );
 
erp.trial = erp.trial(:,:,nearest(erp.time,0):nearest(erp.time,1.1));
lv_pretty_errorbar(erp.time(nearest(erp.time,0):nearest(erp.time,1.1)), squeeze(erp.trial), ...
    (squeeze(erp.trial).*0)+0.5, 1);

%% regression of features to find the transformation from one to another
% zscoring does that so they are all on the same scaling..
% 2d data
mu = [0 0];
Sigma = [1 1.5; 1.5 3]; rng(1)
trn_c1 = mvnrnd(mu,Sigma,100); %outlier = find(sum(trn_c1,2)==max(sum(trn_c1,2),[],1)); trn_c1(outlier,:)=trn_c1(outlier,:)+[10 90];
trn_c1 = zscore(trn_c1,[],1); 
%figure, plot(trn_c1(:,1),trn_c1(:,2),'+')
new_feat2 = (trn_c1(:,1)\trn_c1(:,2))*trn_c1(:,1); % will use this instead of the second feature
hold on, plot(trn_c1(:,1), new_feat2)

% compare the feat1 feat2 relationship to feat1 new_feat2 relationship
figure, plot(trn_c1(:,1),'+')
hold on, plot(trn_c1(:,2),'o')
figure, plot(trn_c1(:,1),'+') % now they are aligned and the value of new_feat2 is approximation of the actual feat2 value ..
% this could help align a feature from different datasets to get the
% abstract relationship between datasets and reduce impact of scaling..
hold on, plot(new_feat2,'o')
% sort the feture values so that the trial to trial calculation is somewhat
% fair and then reduce the datasets to the min of the two and then align
% the feature based on that abstract relationship...
% now between datasets:
% 2d data
mu = [0 0];
Sigma = [1 1.5; 1.5 3]; rng(1)
trn_c1 = mvnrnd(mu,Sigma,100); %outlier = find(sum(trn_c1,2)==max(sum(trn_c1,2),[],1)); trn_c1(outlier,:)=trn_c1(outlier,:)+[10 90];
trn_c2 = trn_c1-5;

trn = [trn_c1;trn_c2];
trn = zscore(trn,[],1);
trn_c1 = trn(1:100,:); trn_c2 = trn(101:end,:);

figure, plot(trn_c1(:,1),trn_c1(:,2),'+')
hold on, plot(trn_c2(:,1),trn_c2(:,2),'o')

rotation = -60; % change rotation to test new combinations .. should be negative for  the features to take the same side for ddifferent classes 
%o.w. the effect means that the classes are swapped and there is no classless rotation that will solve that and this is a problem with the data itself
T = [cosd(rotation) -sind(rotation);sind(rotation) cosd(rotation)];
tst = T*trn';
tst_c1 = tst(:,1:100)'; tst_c2 = tst(:,101:end)';

tst = [tst_c1;tst_c2];
% tst = tst(randperm(size(tst,1)),:); % randomly shuffling
tst = zscore(tst,[],1); 

tst_c1 = tst(1:100,:); tst_c2 = tst(101:end,:);
plot(tst_c1(:,1),tst_c1(:,2),'+')
hold on, plot(tst_c2(:,1),tst_c2(:,2),'o')

% align features:
trn1 = sort(trn(:,1)); tst1 = sort(tst(:,1));
w1 = tst1\trn1;
figure, plot(trn_c1(:,1),trn_c1(:,2),'+')
hold on, plot(trn_c2(:,1),trn_c2(:,2),'o')
trn2 = sort(trn(:,2)); tst2 = sort(tst(:,2));
w2 = tst2\trn2;
plot(tst_c1(:,1)*w1,tst_c1(:,2)*w2,'+')
hold on, plot(tst_c2(:,1)*w1,tst_c2(:,2)*w2,'o')

%% EEG to CNN topos and extracting features with CNN
% from eeg to images of topos for cnn
% this takes trls_image_time .. and because images will be 2dim itself
% this will make a huge array so we need to work on every time pt alone
% so we will loop on time and prepare for cnn and classify like that on
% every time pt and aggregate the result
cfg=[]; cfg.latency=[0.8 1.1];
im = ft_selectdata(cfg, im);

cfg = [];
cfg.layout = 'easycapM1.mat';
lay = ft_prepare_layout(cfg);
lay.label = lower(lay.label); labels = lower(im.label);  % because the labels were sometimes capital and sometimes small
for i=1:length(labels)
    id = ismember(lay.label,labels(i)); temp_lay.label(i)=lay.label(id); temp_lay.pos(i,:)=lay.pos(id,:); temp_lay.height(i)=lay.height(id);
    temp_lay.width(i)=lay.width(id);
end
temp_lay.outline=lay.outline; temp_lay.mask=lay.mask;
for i=1:size(im.trial,3)
    cfg=[];
    cfg.data=im;
    cfg.data.trial = cfg.data.trial(:,:,i); % working on the current time pt
    cfg.method = 'eeg_to_topos_video';
    cfg.data.layout = temp_lay;
    topos_data = lv_feature_extractor(cfg);

    % replacing edge nans with a baseline which is the mean of the image
    baseline=nanmean(topos_data.trial,2);
    topos_data.trial(:, isnan(mean(topos_data.trial,1))) = repmat(baseline,1, sum(isnan(mean(topos_data.trial,1))) );
    % we have the images in the 2nd dim now we use cnn for classification
    cfg=[];
    cfg.classifier_type = {'CNN'};
    cfg.do_parallel = 0; cfg.method = 'timextime';
    cfg.perf_measure = 'auc';%'dval'; % returns trnTime_tstTrls_tstTime .. because every trn time gives the whole tst set(trlsxtime)
    cfg.trn.trial(:,:,1) = topos_data.trial;
    cfg.folds=5;
    cfg.trn.trialinfo=clabel_hand_train;
    txtAuc(nn,:,i) = lv_classify(cfg);
end
% assuming  that we checked the curve and the max was at time 0.8
cfg=[]; cfg.latency=0.8;
im_pt = ft_selectdata(cfg, im);

cfg=[]; % converting peak time pt to topo.
cfg.data=im_pt;
cfg.method = 'eeg_to_topos_video';
cfg.data.layout = temp_lay;
peak_topo = lv_feature_extractor(cfg);
% replacing nans with a baseline which is the mean of the image
baseline=nanmean(peak_topo.trial,2);
peak_topo.trial(:, isnan(mean(peak_topo.trial,1))) = repmat(baseline,1, sum(isnan(mean(peak_topo.trial,1))) );

cfg=[]; % getting the model of the peak
cfg.data.trial=peak_topo.trial;
cfg.data.trialinfo = clabel_hand_train;
cfg.transform = 0;
cfg.method = 'CNN'; cfg.fs=200;
cnn_model = lv_feature_extractor(cfg);
% now we use the data to transform using the learned peak_topo cnn
% so this should be again a loop on all topos of all trials and time
% pts but now they will be reduced using cnn feature extraction
data_after_cnn=[];
for i=1:size(im.trial,3)
    cfg=[];
    cfg.data=im;
    cfg.data.trial = cfg.data.trial(:,:,i); % working on the current time pt
    cfg.method = 'eeg_to_topos_video';
    cfg.data.layout = temp_lay;
    topos_data = lv_feature_extractor(cfg);

    % replacing nans with a baseline which is the mean of the image
    baseline=nanmean(topos_data.trial,2);
    topos_data.trial(:, isnan(mean(topos_data.trial,1))) = repmat(baseline,1, sum(isnan(mean(topos_data.trial,1))) );

    cfg=[]; % transforming with the model of the peak
    cfg.data.trial=peak_topo.trial;
    cfg.data.trialinfo = clabel_hand_train;
    cfg.data_to_transform.trial = topos_data.trial;
    cfg.transform = 1; cfg.net = cnn_model;
    cfg.method = 'CNN'; cfg.fs=200;
    temp = lv_feature_extractor(cfg);
    data_after_cnn.trial(:,:,i) = temp.trial;
end
% if we have two datasets like sleep and wake we will need to transform
% both and then do the classification
% RNN or another classifier
cfg=[];
cfg.classifier_type = {'lda'};
cfg.do_parallel = 0; cfg.method = 'timextime';
cfg.perf_measure = 'auc';%'dval'; % returns trnTime_tstTrls_tstTime .. because every trn time gives the whole tst set(trlsxtime)
cfg.trn.trial(:,:,1) = data_after_cnn.trial;
cfg.folds=5;
cfg.trn.trialinfo=clabel_hand_train;
txtAuc = lv_classify(cfg);

cfg=[];
cfg.method = 'deep';
cfg.classifier_type = {'RNN'};
cfg.do_parallel = 0;
cfg.trn.trial = cat(3,data_after_cnn.trial, data_after_cnn.trial); % put the actual data here
cfg.trn.trialinfo=clabel_hand_train;  cfg.folds=5;
aucrnn = lv_classify(cfg);

% PCA spatial patterns and then CNN
pcs = cfg.eigVects(cfg.chosen);
for i=1:sum(cfg.chosen)
    for j=1:size(result_tst.trial,1)
        temp(j,:) = pcs(:,i)' * cov(squeeze(result_tst.trial(j,:,:)))'; % trl_ch
    end
    % CNN for this component as if it is the time pt we want
    % then put all the transformation of the datasets here but don't put
    % the classification with RNN/lda this should be outside when we
    % aggregate the transformed data..

end
 


