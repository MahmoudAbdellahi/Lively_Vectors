% opens scripts from a common path in a loop and hard coded paths if
% special scripts

%% special
edit('D:\sul''s code\Matt\sleep\erps\Organised\start_up.m')
edit('D:\sul''s code\Matt\sleep\preproc\segment_imagery_data.m')
edit('D:\sul''s code\Matt\sleep\preproc\segment_sleep_data.m')
edit('D:\sul''s code\Matt\sleep\erps\Organised\Miguel\measures\fn_closeloop_compute_SOSpindle_multiplemeasures.m')
%% loop
% paths 
newexp_path = 'D:\sul''s code\Matt\sleep\erps\Organised\';
% file names
names_list = {...
    'REM_ERP_txt_sleepWake',...
    'REM_ERP_txt_sleepWake_aligned',...
    'feat_measure_classifier'...
    'feat_measure_classifier_CorrectVsIncorrect'...    
    'cal_wake_toi'...
    'hilbert_On_correct'...
    'PHASE_analysis_hilbert'... 
    'align_erps'
    }; 
for i=1:length(names_list), edit([newexp_path names_list{i} '.m']); end

%% loop
% paths 
newexp_path = 'D:\sul''s code\Matt\sleep\erps\Organised\New exp\';
% file names
names_list = {...
    'hilbert_On_correct_newexp',...
    'PHASE_analysis_hilbert_newexp'...
    'fn_closeloop_compute_SOSpindle_multiplemeasures_newexp'...
    'feat_measure_classifier_CorrectVsIncorrect_newexp'...
    'ERP_Analysis_behavior_corrSWS_newexp'...
    'addDirs',...
    'lv_pipeline'...
    'lv_clean_segmented'...
    'lv_segment_filter_raw'...
    'lv_segment_filter_raw_lateral'...
    'lv_manual_cleaning'...
    'lv_plot_topo'...
    'lv_feature_extractor'...
    'lv_erp'...
    'lv_tf'...
    'lv_detect_phenomenon'...
    'lv_classify'...
    'lv_post_classification'...
    'lv_reduce_trials'...  
    'lv_build_pipeline'...
    'lv_component_analysis'...
    'lv_examples'...
    'lv_align'... 
    'lv_slider'... 
    'lv_fix_channels'... 
    'lv_align_pipeline'... 
    'lv_trlxtrl_similarity'... 
    'lv_dim_handle'... 
    'lv_EEG_to_CNN'... 
    }; 
for i=1:length(names_list), edit([newexp_path names_list{i} '.m']); end



 








