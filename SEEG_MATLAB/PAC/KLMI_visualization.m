%% ============================================================
% KLMI PAC visualization
%
% Description:
%   Generate KLMI comodulograms for FLC and hippocampus
%   across different experimental states.
%
% Input:
%   KLMI_CONTINUOUS_*.mat
%
% Processing:
%   - Subject averaging
%   - Gaussian smoothing
%   - KLMI visualization
%
% MATLAB version:
%   R2020b
%
% ============================================================


clc;
clear;
close all;


%% ========= Path =========

% Folder containing KLMI_CONTINUOUS_*.mat

dataDir = 'YOUR_KLMI_RESULT_FOLDER';


figDir = fullfile(dataDir,'figures');


if ~exist(figDir,'dir')
    mkdir(figDir);
end



%% ========= State =========

stateNames  = {'MAIN','ANE','ETT','AWA','WARD'};

stateLabels = {'ANE','R-5','R+5','R+25','WARD'};

stateFiles  = stateNames;



%% ========= ROI =========

roi_names = {'FLC','HPC'};



%% ========= Smoothing parameter =========

sigma = 1.2;



%% ========= Main loop =========

for s = 1:length(stateNames)


    state = stateNames{s};


    fprintf('Processing %s...\n',state);



    %% ===== Load =====

    file = fullfile(dataDir,...
        ['KLMI_CONTINUOUS_' state '.mat']);



    if ~exist(file,'file')

        warning('Missing file: %s',file);

        continue;

    end



    load(file);   
    % contains:
    % klmi_FLC
    % klmi_HPC
    % lf_lower
    % hf_lower



    %% ===== ROI loop =====

    for r = 1:length(roi_names)


        if r == 1

            data_cell = klmi_FLC;

        else

            data_cell = klmi_HPC;

        end



        %% ===== Concatenate subjects =====

        valid_idx = ~cellfun(@isempty,data_cell);



        if sum(valid_idx)==0

            warning('%s - %s empty',...
                state,...
                roi_names{r});

            continue;

        end



        tmp = cat(3,data_cell{valid_idx});



        %% ===== Average =====

        img = mean(tmp,3);



        %% ===== Gaussian smoothing =====

        img_smooth = imgaussfilt(img,sigma);



        %% ================= Plot =================


        figure('Color','w',...
            'Position',[200 200 650 520]);



        imagesc(lf_lower,...
                hf_lower,...
                img_smooth');

        axis xy



        caxis([0 8e-4]);

        colormap(turbo)



        cb=colorbar;

        cb.FontSize=12;



        set(gca,...
            'FontSize',14,...
            'LineWidth',1.2);



        xlabel('Low Frequency (Hz)',...
            'FontWeight','bold');


        ylabel('High Frequency (Hz)',...
            'FontWeight','bold');



        title([stateLabels{s} ' - ' roi_names{r}],...
            'FontSize',16,...
            'FontWeight','bold');



        %% ===== Save =====

        exportgraphics(gcf,...
            fullfile(figDir,...
            ['KLMI_' stateFiles{s} '_' ...
            roi_names{r} '.tif']),...
            'Resolution',300);



        close

    end

end



disp('===== ALL FIGURES DONE =====');