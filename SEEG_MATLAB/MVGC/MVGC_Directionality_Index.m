%% ======================== DI Spectrum (LOG SCALE) ========================

clear; clc; close all;


%% ===== 全局字体 =====

set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultLegendFontName','Times New Roman');


figure('Color','w','Position',[100 100 900 600]); 
hold on;


plotHandle = gobjects(length(state_names),1);



for s = 1:length(state_names)


    thisDir = dirs.(state_names{s});

    files = dir(fullfile(thisDir,'MVGC_*.mat'));


    allDI = [];



    for fIdx = 1:length(files)


        try

            load(fullfile(thisDir,files(fIdx).name),'result');


            GC = result.GC;
            freq = result.freq;



            %% ===== 频率限制 =====

            freq_idx = freq >= 1 & freq <= 100;

            freq_plot = freq(freq_idx);



            %% ===== 两个方向 =====
            %
            % GC矩阵:
            % 1 = FLC
            % 2 = HPC
            %
            % GC(1,2): FLC → HPC
            % GC(2,1): HPC → FLC


            gc_F2H = squeeze(GC(1,2,:)); % FLC → HPC

            gc_H2F = squeeze(GC(2,1,:)); % HPC → FLC



            gc_F2H = gc_F2H(freq_idx);

            gc_H2F = gc_H2F(freq_idx);



            %% ===== 防止除0 =====

            denom = gc_H2F + gc_F2H + eps;



            %% ===== Directionality Index =====
            %
            % DI > 0 : HPC → FLC占优势
            % DI < 0 : FLC → HPC占优势


            DI = (gc_H2F - gc_F2H) ./ denom;


            DI = reshape(DI,1,[]);



            allDI = [allDI; DI];



        catch

            continue;

        end


    end



    %% ===== 无数据跳过 =====

    if isempty(allDI)

        warning('%s 无数据', state_names{s});

        continue;

    end




    %% ===== 统计 =====

    meanDI = mean(allDI,1);

    stdDI  = std(allDI,0,1);

    nSub   = size(allDI,1);



    ci95 = 1.96 * stdDI ./ sqrt(nSub);



    upper = meanDI + ci95;

    lower = meanDI - ci95;




    %% ===== 平滑 =====

    meanDI = smoothdata(meanDI,'gaussian',5);

    upper  = smoothdata(upper,'gaussian',5);

    lower  = smoothdata(lower,'gaussian',5);



    x = freq_plot;



    %% ===== 对齐长度 =====

    minLen = min([length(x),...
                  length(upper),...
                  length(lower),...
                  length(meanDI)]);



    x      = x(1:minLen);

    upper  = upper(1:minLen);

    lower  = lower(1:minLen);

    meanDI = meanDI(1:minLen);




    %% ===== CI Ribbon =====

    patch([x fliplr(x)],...
          [upper fliplr(lower)],...
          color.(state_names{s}),...
          'FaceAlpha',0.2,...
          'EdgeColor','none',...
          'HandleVisibility','off');




    %% ===== 均值曲线 =====

    plotHandle(s) = plot(x,...
                         meanDI,...
                         'Color',...
                         color.(state_names{s}),...
                         'LineWidth',...
                         2.5);



end




%% ======================== 图形设置 ========================


xlabel('Frequency (Hz)',...
    'FontSize',16,...
    'FontName','Times New Roman');



ylabel('Directionality Index (DI)',...
    'FontSize',16,...
    'FontName','Times New Roman');



title('HPC ↔ FLC Directionality Index',...
    'FontWeight','bold',...
    'FontSize',20,...
    'FontName','Times New Roman');



xlim([1 100]);

ylim([-1 1]);



set(gca,...
    'XScale','log',...
    'LineWidth',1.5,...
    'TickDir','out',...
    'FontSize',14,...
    'FontName','Times New Roman');



ax = gca;


ax.XTick = [1 2 4 8 13 30 60 100];

ax.XTickLabel = {'1','2','4','8','13','30','60','100'};




%% ===== 关键参考线 =====

yline(0,'k--','LineWidth',1.2);




%% ===== 外观 =====

grid off;

box off;




%% ===== legend =====

validHandles = plotHandle(isgraphics(plotHandle));


lgd = legend(validHandles,...
             state_labels,...
             'Location','eastoutside');


set(lgd,...
    'Box','off',...
    'FontName','Times New Roman',...
    'FontSize',13);