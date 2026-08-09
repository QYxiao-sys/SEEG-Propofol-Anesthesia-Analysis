%% =====================================================
% FLC-HPC Spectral Granger Causality Visualization
%
% This script plots bidirectional spectral GC:
% FLC → HPC
% HPC → FLC
%
% Data format:
% data/
%   MAIN/
%   ANE/
%   ETT/
%   AWA/
%   WARD/
%
% Each folder contains:
% MVGC_*.mat
%
% Each MAT file should contain:
% result.GC
% result.freq
%
%% =====================================================


clear; clc; close all;


%% ===== Font =====

set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultLegendFontName','Times New Roman');

set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12);



%% =====================================================
% Relative data path
%% =====================================================

codeDir = fileparts(mfilename('fullpath'));

dataRoot = fullfile(codeDir,'..','data');



%% ======================== State ========================

state_names  = {'MAIN','ANE','ETT','AWA','WARD'};

state_labels = {'ANE','R-10','R+5','R+25','WARD'};



%% ======================== Color ========================

% FLC → HPC
color_fwd = [0.85 0.33 0.10];

% HPC → FLC
color_rev = [0 0 0.7];



%% =====================================================
% Main loop
%% =====================================================


for s = 1:length(state_names)


    state = state_names{s};


    thisDir = fullfile(dataRoot,state);


    files = dir(fullfile(thisDir,'MVGC_*.mat'));



    allGC_fwd = [];

    allGC_rev = [];



    for fIdx = 1:length(files)


        try

            load(fullfile(thisDir,...
                files(fIdx).name),...
                'result');



            GC = result.GC;

            freq = result.freq;



            %% Frequency range

            freq_idx = freq>=1 & freq<=100;

            freq_plot = freq(freq_idx);



            %% ==============================
            % Direction
            %
            % Channel order:
            % 1 = FLC
            % 2 = HPC
            %
            % GC(1,2):
            % FLC → HPC
            %
            % GC(2,1):
            % HPC → FLC
            %% ==============================


            gc_fwd = squeeze(GC(1,2,:));

            gc_rev = squeeze(GC(2,1,:));



            gc_fwd = gc_fwd(freq_idx);

            gc_rev = gc_rev(freq_idx);



            allGC_fwd(end+1,:) = gc_fwd(:)';

            allGC_rev(end+1,:) = gc_rev(:)';



        catch

            warning('Failed: %s',...
                files(fIdx).name);

            continue

        end

    end



    if isempty(allGC_fwd)

        warning('%s no data',state);

        continue

    end



    %% ================= Statistics =================


    mean_fwd = mean(allGC_fwd,1);

    mean_rev = mean(allGC_rev,1);



    std_fwd = std(allGC_fwd,[],1);

    std_rev = std(allGC_rev,[],1);



    nSub=size(allGC_fwd,1);



    ci_fwd = 1.96*std_fwd./sqrt(nSub);

    ci_rev = 1.96*std_rev./sqrt(nSub);



    upper_fwd = mean_fwd+ci_fwd;

    lower_fwd = mean_fwd-ci_fwd;


    upper_rev = mean_rev+ci_rev;

    lower_rev = mean_rev-ci_rev;



    x=freq_plot;



    %% ================= Plot =================


    figure('Color','w',...
        'Position',[100 100 850 550]);

    hold on;



    %% CI

    patch([x fliplr(x)],...
        [upper_fwd fliplr(lower_fwd)],...
        color_fwd,...
        'FaceAlpha',0.2,...
        'EdgeColor','none',...
        'HandleVisibility','off');



    patch([x fliplr(x)],...
        [upper_rev fliplr(lower_rev)],...
        color_rev,...
        'FaceAlpha',0.2,...
        'EdgeColor','none',...
        'HandleVisibility','off');



    %% Curves

    plot(x,mean_fwd,...
        'Color',color_fwd,...
        'LineWidth',2,...
        'DisplayName',...
        'FLC → HPC');


    plot(x,mean_rev,...
        'Color',color_rev,...
        'LineWidth',2,...
        'DisplayName',...
        'HPC → FLC');



    %% Axis

    xlabel('Frequency (Hz)');

    ylabel('Spectral Granger Causality');


    title(state_labels{s},...
        'FontWeight','bold');



    xlim([1 100]);

    ylim([0 0.12]);


    set(gca,...
        'XScale','log',...
        'LineWidth',1.2,...
        'TickDir','out');



    ax=gca;

    ax.XTick=[1 2 4 8 13 30 60 100];

    ax.XTickLabel={'1','2','4','8','13','30','60','100'};



    box off;

    grid off;



    legend('Location','northeast');



end