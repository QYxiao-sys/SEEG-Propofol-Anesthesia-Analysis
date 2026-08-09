%% ============================================================
% Continuous SEEG KLMI phase-amplitude coupling analysis
%
% Description:
%   Calculate channel-level KLMI PAC during different
%   experimental states and obtain ROI-level PAC values.
%
% States:
%   MAIN, ANE, ETT, AWA, WARD
%
% ROIs:
%   PFC
%   HPC
%
% MATLAB version:
%   R2020b
%
% Required toolbox:
%   EEGLAB
%
% Required function:
%   get_klmi_jia.m
%
% ============================================================


clc;
clear;
close all;



%% ========= 0. Output directory =========


% User-defined output folder

outDir = 'YOUR_OUTPUT_FOLDER';


if ~exist(outDir,'dir')

    mkdir(outDir);

end



%% ========= 1. States =========


stateNames = {'MAIN','ANE','ETT','AWA','WARD'};


% User-defined data root folder

baseDir = 'YOUR_SEEG_DATA_ROOT';



%% ========= 2. Frequency parameters =========


lf_lower  = 0.5:0.25:13;

lf_higher = lf_lower + 0.5;


hf_lower  = 15:2.5:100;

hf_higher = hf_lower + 6;


nbins = 18;



%% =====================================================
%% ================= State loop ========================
%% =====================================================


for sIdx = 1:length(stateNames)


    state = stateNames{sIdx};


    fprintf('\n===== Processing STATE: %s =====\n',...
        state);



    dataDir = fullfile(baseDir,state);



    if ~exist(dataDir,'dir')

        warning('Directory does not exist: %s',...
            dataDir);

        continue;

    end



    files = dir(fullfile(dataDir,'*.set'));


    fileNames = {files.name};



    if isempty(fileNames)

        warning('%s contains no data',state);

        continue;

    end



    %% ===== Subject-level results =====


    klmi_PFC = cell(length(fileNames),1);

    klmi_HPC = cell(length(fileNames),1);



    %% =====================================================
    %% ================= Subject loop ======================
    %% =====================================================


    for m = 1:length(fileNames)


        fprintf('\nSubject %d / %d : %s\n',...
            m,...
            length(fileNames),...
            fileNames{m});



        EEG = pop_loadset(...
            'filename',fileNames{m},...
            'filepath',dataDir);



        EEG = eeg_checkset(EEG);



        %% ===== Continuous data check =====


        if ndims(EEG.data) ~= 2

            error('Dataset is not continuous: %s',...
                fileNames{m});

        end



        chanNames = {EEG.chanlocs.labels};


        nCh = size(EEG.data,1);



        %% ===== ROI index =====


        idx_PFC = find(~cellfun(@isempty,...
            regexp(chanNames,'^PFC')));



        idx_HPC = find(~cellfun(@isempty,...
            regexp(chanNames,'^HPC')));



        %% ===== Initialization =====


        klmi_allCh = zeros(...
            length(lf_lower),...
            length(hf_lower),...
            nCh);



        %% =====================================================
        %% =============== Channel-level PAC ==================
        %% =====================================================


        for ch = 1:nCh


            fprintf('  Channel %d / %d\n',...
                ch,nCh);



            data_ch = double(EEG.data(ch,:));


            data_ch = data_ch(:)';



            %% ===== LF phase =====


            phase_lf = zeros(...
                length(data_ch),...
                length(lf_lower));



            for i = 1:length(lf_lower)


                tmp = eegfilt(...
                    data_ch,...
                    EEG.srate,...
                    lf_lower(i),...
                    lf_higher(i),...
                    0,...
                    [],...
                    0,...
                    'fir1',...
                    0);



                phase_lf(:,i) = ...
                    angle(hilbert(tmp))';


            end



            %% ===== HF amplitude =====


            amplitude_hf = zeros(...
                length(data_ch),...
                length(hf_lower));



            for j = 1:length(hf_lower)


                tmp = eegfilt(...
                    data_ch,...
                    EEG.srate,...
                    hf_lower(j),...
                    hf_higher(j),...
                    0,...
                    [],...
                    0,...
                    'fir1',...
                    0);



                amplitude_hf(:,j)=...
                    abs(hilbert(tmp))';


            end



            %% ===== KLMI =====


            out = get_klmi_jia(...
                phase_lf,...
                amplitude_hf,...
                nbins);



            tmpMI = out.MI;



            if ndims(tmpMI)==3

                tmpMI=squeeze(tmpMI);

            end



            klmi_allCh(:,:,ch)=tmpMI;



        end



        %% =====================================================
        %% ================= ROI average ======================
        %% =====================================================



        % ===== PFC =====

        if ~isempty(idx_PFC)


            klmi_PFC{m}=...
                mean(klmi_allCh(:,:,idx_PFC),3);


        else


            warning('No PFC in %s',...
                fileNames{m});


            klmi_PFC{m}=[];


        end



        % ===== HPC =====

        if ~isempty(idx_HPC)


            klmi_HPC{m}=...
                mean(klmi_allCh(:,:,idx_HPC),3);


        else


            warning('No HPC in %s',...
                fileNames{m});


            klmi_HPC{m}=[];


        end



        clear EEG klmi_allCh


    end



    %% =====================================================
    %% ================= Save ===============================
    %% =====================================================


    save(fullfile(outDir,...
        ['KLMI_CONTINUOUS_' state '.mat']),...
        'klmi_PFC',...
        'klmi_HPC',...
        'lf_lower',...
        'hf_lower',...
        '-v7.3');



    fprintf('===== %s DONE =====\n',state);


end



disp('===== ALL STATES FINISHED =====');