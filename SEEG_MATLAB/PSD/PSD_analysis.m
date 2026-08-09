%% ============================================================
% Subject-level PSD analysis for SEEG recordings
%
% Description:
%   Compute power spectral density (PSD) for predefined ROIs
%   across different experimental states.
%
% MATLAB version:
%   R2020b
%
% Required toolbox:
%   EEGLAB
%
% ============================================================


clear; clc; close all;


%% ===== 加载 EEGLAB =====

% User-defined EEGLAB path
eeglab_path = 'YOUR_EEGLAB_PATH';

addpath(genpath(eeglab_path));

eeglab;   % 初始化（必须）


%% ====================================================
%% ================= 参数 ==============================
%% ====================================================

Fs   = 512;

L    = 2048;

NFFT = 2^nextpow2(L);


f = Fs/2 * linspace(0,1,NFFT/2+1);

freq_idx  = f >= 0.5 & f <= 100;

freq_plot = f(freq_idx);



%% ====================================================
%% ================= 状态 ==============================
%% ====================================================


state_map.names  = {'MAIN','ANE','ETT','AWA','WARD'};

state_map.labels = {'ANE','R-5','R+5','R+25','WARD'};


% User-defined data root path
data_root = 'YOUR_DATA_ROOT';


state_map.paths = { ...
    fullfile(data_root,'MAIN'), ...
    fullfile(data_root,'ANE'), ...
    fullfile(data_root,'ETT'), ...
    fullfile(data_root,'AWA'), ...
    fullfile(data_root,'WARD')};


states = state_map.names;



%% ====================================================
%% ================= ROI ===============================
%% ====================================================

ROI_list = {'PFC','HPC'};



%% ====================================================
%% ================= 输出路径 ===========================
%% ====================================================


save_dir = 'YOUR_OUTPUT_FOLDER';


if ~exist(save_dir,'dir')

    mkdir(save_dir);

end



%% ====================================================
%% ================= 初始化 =============================
%% ====================================================


for r = 1:length(ROI_list)

    PSD_ROI.(ROI_list{r}) = struct();

end



disp('======================================')

disp('      SUBJECT LEVEL PSD COMPUTE')

disp('======================================')



%% ====================================================
%% ================= 主循环 ============================
%% ====================================================


for s = 1:length(states)


    fprintf('\n===== STATE: %s =====\n',states{s});


    files = dir(fullfile(state_map.paths{s}, '*.set'));


    ROI_subject.PFC = [];

    ROI_subject.HPC = [];



    %% ---------- 每个被试 ----------

    for subj = 1:length(files)


        fprintf('Subject %d / %d\n',...
            subj,length(files));


        EEG = pop_loadset(...
            'filename',files(subj).name,...
            'filepath',state_map.paths{s});


        labels = {EEG.chanlocs.labels};



        %% ===== ROI indexing =====


        idx.PFC = find(strncmpi(labels,'PFC',3));


        % HPC = TH + Tb

        idx.HPC = find(~cellfun('isempty',strfind(labels,'HPC_TH')) | ...
                       ~cellfun('isempty',strfind(labels,'HPC_Tb')));



        nEpoch = EEG.trials;


        subj_data.PFC = [];

        subj_data.HPC = [];


        ROI_names = fieldnames(idx);



        %% ---------- PSD ----------

        for r = 1:length(ROI_names)


            roi_name = ROI_names{r};

            roi_idx  = idx.(roi_name);



            if isempty(roi_idx)

                continue;

            end



            for ch = roi_idx


                for ep = 1:nEpoch


                    y = double(squeeze(EEG.data(ch,:,ep)));


                    Pxx = pwelch(y,...
                                 hamming(L),...
                                 L/2,...
                                 NFFT,...
                                 Fs);



                    subj_data.(roi_name)(end+1,:) = ...
                        10*log10(Pxx(freq_idx)' + eps);



                end

            end

        end



        %% ===== SUBJECT 平均 =====

        for r = 1:length(ROI_names)


            roi_name = ROI_names{r};



            if isempty(subj_data.(roi_name))

                continue;

            end



            ROI_subject.(roi_name)(end+1,:) = ...
                mean(subj_data.(roi_name),1);


        end


    end



    %% ===== 保存该状态 =====

    for r = 1:length(ROI_list)


        roi_name = ROI_list{r};


        PSD_ROI.(roi_name).(states{s}) = ...
            ROI_subject.(roi_name);


    end


end



disp('===== SUBJECT PSD COMPUTATION FINISHED =====')



%% ====================================================
%% ================= 保存 MAT 文件 =====================
%% ====================================================


disp('===== SAVING ROI MAT FILES =====')



for r = 1:length(ROI_list)


    roi_name = ROI_list{r};


    PSD   = PSD_ROI.(roi_name);

    freq  = freq_plot;



    save_name = fullfile(save_dir,...
        ['PSD_' roi_name '.mat']);



    % 保存 labels

    state_labels = state_map.labels;



    save(save_name,...
        'PSD',...
        'freq',...
        'states',...
        'state_labels',...
        '-v7.3');



    fprintf('Saved -> %s\n',save_name);



end



disp('======================================')

disp('        ALL MAT FILES SAVED')

disp('======================================')