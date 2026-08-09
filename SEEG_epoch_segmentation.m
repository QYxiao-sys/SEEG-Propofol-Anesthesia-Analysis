%% ============================================================
% SEEG epoch segmentation
%
% Description:
%   Segment continuous SEEG recordings into overlapping epochs
%   using EEGLAB.
%
% Current settings:
%   Epoch length: 4 s
%   Step size:    2 s
%   Overlap:      50%
%
% MATLAB version:
%   R2020b
%
% Required toolbox:
%   EEGLAB
%
% ============================================================


clearvars;
clc;


%% ============================
% User-defined paths
% ============================

% Input folder containing continuous SEEG datasets
input_dir = 'YOUR_INPUT_FOLDER';


% Output folder for epoched datasets
output_dir = 'YOUR_OUTPUT_FOLDER';


if ~exist(output_dir,'dir')

    mkdir(output_dir);

end



%% ============================
% Epoch parameters
% ============================

epoch_length = 4;     % seconds

epoch_step = 2;       % seconds

overlap_ratio = ...
    (epoch_length-epoch_step) / epoch_length;



fprintf('Epoch length: %.1f s\n',epoch_length);

fprintf('Step size: %.1f s\n',epoch_step);

fprintf('Overlap: %.1f %%\n',...
    overlap_ratio*100);



%% ============================
% Find datasets
% ============================

files = dir(fullfile(input_dir,'*.set'));


if isempty(files)

    error('No EEGLAB datasets found in %s',input_dir);

end



%% ============================
% Batch epoching
% ============================

for subj = 1:length(files)


    fprintf('\nProcessing %d/%d: %s\n',...
        subj,...
        length(files),...
        files(subj).name);



    %% Load continuous SEEG

    EEG = pop_loadset(...
        'filename',files(subj).name,...
        'filepath',input_dir);



    %% Create overlapping epochs

    EEG = eeg_regepochs(...
        EEG,...
        'recurrence',epoch_step,...
        'limits',[0 epoch_length],...
        'rmbase',NaN);



    EEG = eeg_checkset(EEG);



    %% Save

    output_name = ...
        ['epoch_' files(subj).name];


    pop_saveset(...
        EEG,...
        'filename',output_name,...
        'filepath',output_dir);



    fprintf('Saved: %s\n',output_name);


end



disp('===== SEEG epoch segmentation completed =====');
