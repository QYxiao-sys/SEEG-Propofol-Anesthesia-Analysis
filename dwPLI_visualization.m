%% ============================================================
% dwPLI visualization across experimental states
%
% Description:
%   Visualize FLC-HPC dwPLI spectra across five experimental
%   states with 95% confidence intervals.
%
% Input:
%   dwPLI_STATE.mat
%
% Output:
%   dwPLI spectrum plots
%
% ROI:
%   FLC - Frontal lobe cortex
%   HPC - Hippocampus
%
% MATLAB:
%   R2020b
%
% ============================================================


clear;
clc;
close all;



%% ================= Global font =================


set(0,'DefaultAxesFontName','Times New Roman');

set(0,'DefaultTextFontName','Times New Roman');

set(0,'DefaultLegendFontName','Times New Roman');



%% ================= Path =================


% ===== User configuration =====

dataDir = 'YOUR_dwPLI_RESULT_FOLDER';



%% ================= State order =================


states = {'MAIN',...
          'ANE',...
          'ETT',...
          'AWA',...
          'WARD'};



% Labels used in manuscript

stateLabels = {'ANE',...
               'R-5',...
               'R+5',...
               'R+25',...
               'WARD'};



%% ================= Colors =================


colors = struct( ...
    'MAIN',[0.85 0.33 0.10], ...
    'ANE' ,[0.93 0.69 0.13], ...
    'ETT' ,[0.00 0.45 0.74], ...
    'AWA' ,[0.47 0.67 0.19], ...
    'WARD',[0 0 0]);



%% ================= Figure =================


figure('Color','w',...
       'Position',[200 200 900 600]);


hold on;


h_lines = [];

legend_labels = {};



%% ================= Main loop =================


for s = 1:length(states)


    state = states{s};



    %% ===== Load =====


    file = fullfile(dataDir,...
        ['dwPLI_' state '.mat']);



    if ~exist(file,'file')


        warning('%s file missing',state);

        continue;


    end



    load(file);
    
    % Variables:
    % dwpli_all : frequency × subjects
    % freq
    % subj_names



    %% ===== Frequency =====


    if isstruct(freq)

        freqs = freq.freq(:);

    else

        freqs = freq(:);

    end



    %% ===== Check =====


    if size(dwpli_all,1) ~= length(freqs)


        warning('%s frequency mismatch',state);

        continue;


    end



    %% ================= Statistics =================


    medS = median(dwpli_all,2);


    nSubj = size(dwpli_all,2);



    sem = std(dwpli_all,0,2) ./ sqrt(nSubj);


    CI95 = 1.96 * sem;



    CI_low = medS - CI95;

    CI_high = medS + CI95;



    %% ===== Gaussian smoothing =====


    medS = smoothdata(medS,...
        'gaussian',5);


    CI_low = smoothdata(CI_low,...
        'gaussian',5);


    CI_high = smoothdata(CI_high,...
        'gaussian',5);



    %% ================= CI ribbon =================


    fill([freqs; flipud(freqs)],...
         [CI_low; flipud(CI_high)],...
         colors.(state),...
         'FaceAlpha',0.15,...
         'EdgeColor','none');



    %% ================= Main curve =================


    h = plot(freqs,...
             medS,...
             'Color',colors.(state),...
             'LineWidth',2.5);



    h_lines = [h_lines h];


    legend_labels{end+1}=stateLabels{s};



end



%% ================= Axis =================


xlabel('Frequency (Hz)',...
    'FontSize',16,...
    'FontName','Times New Roman');



ylabel('dwPLI',...
    'FontSize',16,...
    'FontName','Times New Roman');



title('FLC–HPC dwPLI (95% CI)',...
    'FontSize',20,...
    'FontWeight','bold');



set(gca,...
    'XScale','log',...
    'XLim',[1 100],...
    'YLim',[0 0.5],...
    'FontSize',14,...
    'FontName','Times New Roman',...
    'LineWidth',1.5,...
    'TickDir','out');



%% ===== Log frequency ticks =====


ax = gca;


ax.XTick = [1 2 4 8 13 30 60 100];


ax.XTickLabel = ...
    {'1','2','4','8','13','30','60','100'};



box off;

grid off;



%% ================= Legend =================


lgd = legend(h_lines,...
             legend_labels,...
             'Location','northeast');



set(lgd,...
    'FontName','Times New Roman',...
    'FontSize',13);



%% ================= Save =================


print(gcf,...
    fullfile(dataDir,...
    'dwPLI_FLC_HPC_all_states.png'),...
    '-dpng',...
    '-r300');



print(gcf,...
    fullfile(dataDir,...
    'dwPLI_FLC_HPC_all_states.eps'),...
    '-depsc2',...
    '-r300');



disp('===== FLC-HPC dwPLI visualization finished =====');