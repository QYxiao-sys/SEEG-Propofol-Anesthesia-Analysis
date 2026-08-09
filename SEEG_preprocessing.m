%% ============================================================
% SEEG preprocessing pipeline
%
% Description:
%   Preprocessing of SEEG recordings exported in EEGLAB format.
%
% Processing steps:
%   1. Load EEGLAB .set files
%   2. Resample to 512 Hz
%   3. Band-pass filtering (0.1-150 Hz)
%   4. Line noise removal (50 Hz and 100 Hz harmonics)
%   5. Keep the first event marker
%   6. Save processed datasets
%
% MATLAB version:
%   R2020b
%
% Required toolbox:
%   EEGLAB
% ============================================================


clearvars;
clc;


%% ============================
% User-defined paths
% ============================

% Folder containing raw EEGLAB files
input_dir = 'YOUR_INPUT_FOLDER';

% Folder for processed files
output_dir = 'YOUR_OUTPUT_FOLDER';


if ~exist(output_dir,'dir')
    mkdir(output_dir);
end



%% ============================
% Preprocessing parameters
% ============================

params.target_fs = 512;

params.bandpass_low  = 0.1;
params.bandpass_high = 150;

params.notch_freq1 = [49 51];
params.notch_freq2 = [99 101];



%% ============================
% Find input files
% ============================

files = dir(fullfile(input_dir,'*.set'));


if isempty(files)

    error('No EEGLAB .set files found in: %s', input_dir);

end



%% ============================
% Batch preprocessing
% ============================

for i = 1:length(files)


    fprintf('\nProcessing %d/%d: %s\n',...
        i,length(files),files(i).name);



    %% Load dataset

    EEG = pop_loadset(...
        'filename',files(i).name,...
        'filepath',input_dir);



    %% Resampling

    EEG = pop_resample(...
        EEG,...
        params.target_fs);



    %% Band-pass filtering

    EEG = pop_eegfiltnew(...
        EEG,...
        params.bandpass_low,...
        params.bandpass_high);



    %% 50 Hz notch filtering

    EEG = pop_eegfiltnew(...
        EEG,...
        params.notch_freq1(1),...
        params.notch_freq1(2),...
        [],...
        1);



    %% 100 Hz harmonic notch filtering

    EEG = pop_eegfiltnew(...
        EEG,...
        params.notch_freq2(1),...
        params.notch_freq2(2),...
        [],...
        1);



    %% Keep first event only

    if isfield(EEG,'event') && length(EEG.event)>1

        EEG.event = EEG.event(1);

        EEG = eeg_checkset(EEG);

    end



    %% Save processed dataset

    EEG = pop_saveset(...
        EEG,...
        'filename',files(i).name,...
        'filepath',output_dir);



    fprintf('Completed: %s\n',files(i).name);


end



disp('All SEEG preprocessing completed.');
