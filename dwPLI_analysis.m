%% ============================================================
% dwPLI functional connectivity analysis using FieldTrip
%
% Description:
%   Compute debiased weighted phase lag index (dwPLI)
%   between frontal lobe cortex (FLC) and hippocampus (HPC)
%   across different experimental states.
%
% Method:
%   - FieldTrip mtmfft
%   - DPSS multitaper
%   - wPLI debiased estimator
%
% Frequency range:
%   1-100 Hz
%
% Output:
%   Subject-level FLC-HPC dwPLI spectra
%
% MATLAB:
%   R2020b
%
% ============================================================


clear;
clc;

try
    close all;
catch
end



%% ============================================================
%% ================= FieldTrip initialization =================
%% ============================================================


% ===== User configuration =====

fieldtrip_path = 'YOUR_FIELDTRIP_PATH';


restoredefaultpath;

addpath(fieldtrip_path);

ft_defaults;



%% ============================================================
%% ================= Basic parameters =========================
%% ============================================================


% Folder containing epoch .set files

data_root = 'YOUR_EPOCH_DATA_FOLDER';



states = {'ANE',...
          'AWA',...
          'ETT',...
          'MAIN',...
          'WARD'};



% Output folder

output_dir = 'YOUR_OUTPUT_FOLDER';



if ~exist(output_dir,'dir')

    mkdir(output_dir);

end



%% ============================================================
%% ================= Multitaper parameters ====================
%% ============================================================


cfg_freq = [];

cfg_freq.method      = 'mtmfft';

cfg_freq.output      = 'powandcsd';

cfg_freq.foilim      = [1 100];

cfg_freq.tapsmofrq   = 2;

cfg_freq.taper       = 'dpss';

cfg_freq.keeptrials  = 'yes';

cfg_freq.pad         = 'nextpow2';



results = struct();



%% ============================================================
%% ================= Main analysis loop =======================
%% ============================================================


for s = 1:length(states)


    state_name = states{s};


    fprintf('\n========== State: %s ==========\n',...
        state_name);



    state_path = fullfile(data_root,state_name);


    set_files = dir(fullfile(state_path,'*.set'));



    subj_results = {};

    subj_names   = {};



    %% ================= Subject loop =================


    for subj = 1:length(set_files)



        filename = fullfile(...
            state_path,...
            set_files(subj).name);



        fprintf('Processing: %s\n',...
            set_files(subj).name);



        %% ===== 1. Load data with FieldTrip =====


        cfg = [];

        cfg.dataset = filename;


        data_ft = ft_preprocessing(cfg);



        %% ===== Trial quality check =====


        if length(data_ft.trial) < 5


            warning('Too few trials, skipped');

            continue;


        end



        %% =====================================================
        %% ===== 2. Select FLC and HPC channels ===============
        %% =====================================================


        labels = data_ft.label;



        idx_FLC = find(...
            startsWith(labels,...
            'FLC',...
            'IgnoreCase',true));



        idx_HPC = find(...
            startsWith(labels,...
            'HPC',...
            'IgnoreCase',true));



        if isempty(idx_FLC) || isempty(idx_HPC)


            warning('Skip: missing FLC or HPC channels');

            continue;


        end



        %% =====================================================
        %% ===== 3. Frequency analysis =========================
        %% =====================================================


        freq = ft_freqanalysis(...
            cfg_freq,...
            data_ft);



        %% =====================================================
        %% ===== 4. Construct FLC-HPC channel pairs ============
        %% =====================================================


        channelcmb = {};



        for f = 1:length(idx_FLC)


            for h = 1:length(idx_HPC)



                channelcmb{end+1,1}=...
                    labels{idx_FLC(f)};



                channelcmb{end,2}=...
                    labels{idx_HPC(h)};



            end


        end



        %% =====================================================
        %% ===== 5. dwPLI calculation ==========================
        %% =====================================================


        cfg_conn = [];

        cfg_conn.method = 'wpli_debiased';

        cfg_conn.channelcmb = channelcmb;



        conn = ft_connectivityanalysis(...
            cfg_conn,...
            freq);



        if isfield(conn,...
                'wpli_debiasedspctrm')


            dwpli_mat = ...
                conn.wpli_debiasedspctrm;


        else


            error('dwPLI output not found');


        end



        %% =====================================================
        %% ===== 6. ROI-level summary ==========================
        %% =====================================================


        % Median across all FLC-HPC channel pairs

        dwpli_subj = median(...
            dwpli_mat,...
            1)';



        %% =====================================================
        %% ===== 7. Store subject result =======================
        %% =====================================================


        subj_results{end+1}=dwpli_subj;


        subj_names{end+1}=...
            set_files(subj).name;



    end



    %% =====================================================
    %% ===== 8. Save state-level result =====================
    %% =====================================================


    if ~isempty(subj_results)



        % frequency × subjects

        dwpli_all = cat(2,...
            subj_results{:});



        results.(state_name).dwpli_all = ...
            dwpli_all;



        results.(state_name).freq = ...
            freq.freq;



        results.(state_name).subjects = ...
            subj_names;



        fprintf('Valid subjects: %d\n',...
            size(dwpli_all,2));



        save(fullfile(output_dir,...
            ['dwPLI_' state_name '.mat']),...
            'dwpli_all',...
            'subj_names',...
            'freq');



    else


        fprintf('No valid subjects\n');


    end



end



%% ============================================================
%% ================= Save all states ==========================
%% ============================================================


save(fullfile(output_dir,...
    'dwPLI_all_states.mat'),...
    'results');



disp('======================================');

disp('FLC-HPC dwPLI computation finished');

disp('======================================');