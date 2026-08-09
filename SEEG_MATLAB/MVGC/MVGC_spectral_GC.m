%% ============================================================
% SEEG ROI spectral Granger causality using MVGC toolbox
%
% Description:
%   Compute frequency-domain Granger causality between
%   frontal lobe cortex (FLC) and hippocampus (HPC)
%
% Input:
%   Epoch EEG (.set)
%
% Output:
%   MVGC_*.mat
%
% ROI:
%   1 - FLC
%   2 - HPC
%
% GC direction:
%   GC(1,2,:) : FLC -> HPC
%   GC(2,1,:) : HPC -> FLC
%
% MATLAB:
%   R2020b
%
% ============================================================


clear;
clc;
close all;



%% ================= Parameters =================


Fs = 256;

nFreq = 128;

maxOrder = 20;



%% ================= User paths =================


dataDir = 'YOUR_EPOCH_DATA_FOLDER';


outDir = 'YOUR_OUTPUT_FOLDER';



%% ================= Toolbox =================


addpath(genpath('YOUR_MVGC_TOOLBOX_FOLDER'));

addpath(genpath('YOUR_EEGLAB_FOLDER'));



if ~exist(outDir,'dir')

    mkdir(outDir);

end



files = dir(fullfile(dataDir,'*.set'));



%% ================= Frequency axis =================


freqVec = linspace(0,...
                   Fs/2,...
                   nFreq);



%% ================= Main loop =================


for fIdx = 1:length(files)


    fprintf('\n========== %s ==========\n',...
        files(fIdx).name);



    try


        %% ===== Load EEG =====


        EEG = pop_loadset(...
            'filename',...
            files(fIdx).name,...
            'filepath',...
            dataDir);



        data = double(EEG.data);


        chanLabels = string({EEG.chanlocs.labels});



        %% ===== ROI selection =====


        FLC_idx = find(...
            startsWith(chanLabels,...
            'FLC',...
            'IgnoreCase',true));



        HPC_idx = find(...
            startsWith(chanLabels,...
            'HPC',...
            'IgnoreCase',true));



        fprintf('FLC: %d   HPC: %d\n',...
            length(FLC_idx),...
            length(HPC_idx));



        %% ===== ROI average =====


        roiData = [];

        roiNames = {};



        if ~isempty(FLC_idx)


            roiData(1,:,:) = ...
                squeeze(mean(data(FLC_idx,:,:),1));


            roiNames{end+1}='FLC';


        end



        if ~isempty(HPC_idx)


            roiData(end+1,:,:) = ...
                squeeze(mean(data(HPC_idx,:,:),1));


            roiNames{end+1}='HPC';


        end



        if size(roiData,1)<2


            warning('ROI number <2, skip');

            continue;


        end



        fprintf('ROI dimension: %d × %d × %d\n',...
            size(roiData));



        %% ================= Preprocessing =================


        for ch=1:size(roiData,1)


            tmp=squeeze(roiData(ch,:,:));


            tmp=detrend(tmp,'constant');


            tmp=zscore(tmp,0,1);


            roiData(ch,:,:)=tmp;


        end



        if any(isnan(roiData(:)))


            warning('NaN detected');

            continue;


        end



        %% ================= VAR order =================


        [~,BIC]=tsdata_to_infocrit(...
            roiData,...
            maxOrder,...
            'OLS');



        [~,order]=min(BIC);


        order=min(order,8);



        fprintf('VAR order = %d\n',...
            order);



        %% ================= VAR model =================


        [A,SIG]=tsdata_to_var(...
            roiData,...
            order,...
            'OLS');



        if isbad(A)


            warning('VAR failed');

            continue;


        end



        SIG=SIG+eye(size(SIG))*1e-6;



        if var_specrad(A)>=1


            warning('VAR unstable');

            continue;


        end



        %% ================= Autocovariance =================


        [G,info]=var_to_autocov(A,SIG);



        if info.error || isbad(G)


            warning('Autocov failed');

            continue;


        end



        %% ================= Spectral GC =================


        fres=nFreq-1;


        GC=autocov_to_spwcgc(G,fres);



        % remove diagonal

        for ch=1:size(GC,1)

            GC(ch,ch,:)=0;

        end



        GC(~isfinite(GC))=0;



        %% ================= Force frequency length =================


        if size(GC,3)~=nFreq


            warning('GC length mismatch');


            GC=GC(:,:,1:nFreq);


        end



        %% ================= Save =================


        result=struct();


        result.GC=GC;


        result.freq=freqVec;


        result.roiNames=roiNames;


        result.subject=files(fIdx).name;



        save(fullfile(outDir,...
            ['MVGC_' files(fIdx).name '.mat']),...
            'result');



        %% ================= QC plot =================


        fig=figure('visible','off');


        plot(freqVec,...
            squeeze(GC(2,1,:)),...
            'r',...
            'LineWidth',1.5);


        hold on;


        plot(freqVec,...
            squeeze(GC(1,2,:)),...
            'k',...
            'LineWidth',1.5);



        xlabel('Frequency (Hz)');

        ylabel('Spectral GC');


        title(files(fIdx).name,...
            'Interpreter','none');



        legend('HPC → FLC',...
               'FLC → HPC');



        grid on;



        saveas(fig,...
            fullfile(outDir,...
            [files(fIdx).name...
            '_FLC_HPC_GC.tif']));



        close(fig);



        fprintf('✔ Finished\n');



    catch ME


        fprintf('✘ Error: %s\n',...
            ME.message);


        continue;


    end


end



fprintf('\n========== MVGC FLC-HPC completed ==========\n');