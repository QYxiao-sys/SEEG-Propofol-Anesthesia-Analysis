%% ============================================================
% dwPLI FLC-HPC band-specific statistics
%
% Description:
%   Band-wise comparison of FLC-HPC dwPLI across anesthesia states
%
% Statistics:
%   Friedman test
%   Paired Wilcoxon signed-rank test
%   Benjamini-Hochberg FDR correction
%   Effect size r
%
% Input:
%   dwPLI_STATE.mat
%
% Output:
%   Band-specific statistical figures
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



%% ================= User parameters =================


data_dir = 'YOUR_dwPLI_RESULT_FOLDER';



states = {'MAIN',...
          'ANE',...
          'ETT',...
          'AWA',...
          'WARD'};



stateLabels = {'ANE',...
               'R-5',...
               'R+5',...
               'R+25',...
               'WARD'};



nStates = length(states);



%% ================= Colors =================


color.MAIN = [0.85 0.33 0.10];

color.ANE  = [0.93 0.69 0.13];

color.ETT  = [0.00 0.45 0.74];

color.AWA  = [0.47 0.67 0.19];

color.WARD = [0 0 0];



state_colors = zeros(nStates,3);


for i = 1:nStates

    state_colors(i,:) = color.(states{i});

end



%% ================= Frequency bands =================


bands = { ...
    'delta',[1 4]; ...
    'theta',[4 8]; ...
    'alpha',[8 13]; ...
    'beta' ,[13 30]; ...
    'gamma',[30 100]};



%% ================= Load data =================


data = struct();


for s = 1:nStates


    tmp = load(fullfile(data_dir,...
        ['dwPLI_' states{s} '.mat']));


    data.(states{s}) = tmp.dwpli_all;


    if isstruct(tmp.freq)

        freq = tmp.freq.freq;

    else

        freq = tmp.freq;

    end


end



%% ================= Common subjects =================


minSubj = inf;


for s = 1:nStates

    minSubj = min(minSubj,...
        size(data.(states{s}),2));

end


fprintf('Common subjects: %d\n',minSubj);



%% ================= Band loop =================


for b = 1:size(bands,1)


    band_name = bands{b,1};

    band_range = bands{b,2};



    fprintf('\n================ %s ================\n',...
        upper(band_name));



    %% ===== Extract frequency =====


    idx = freq >= band_range(1) & ...
          freq <= band_range(2);



    X = zeros(minSubj,nStates);



    for s = 1:nStates


        tmp = data.(states{s});


        tmp = tmp(:,1:minSubj);



        % median dwPLI within frequency band

        band_value = median(tmp(idx,:),1);


        X(:,s)=band_value';



    end



    %% ================= Friedman =================


    p_friedman = friedman(X,1,'off');


    fprintf('Friedman p = %.5f\n',...
        p_friedman);



    %% ================= Pairwise Wilcoxon =================


    pairs = nchoosek(1:nStates,2);


    nPairs=size(pairs,1);



    p_raw=zeros(nPairs,1);

    r_val=zeros(nPairs,1);



    fprintf('\nPairwise comparisons:\n');



    for i=1:nPairs


        s1=pairs(i,1);

        s2=pairs(i,2);



        x1=X(:,s1);

        x2=X(:,s2);



        [p,~,stats]=signrank(x1,x2);



        p_raw(i)=p;



        if isfield(stats,'zval')

            z=stats.zval;

        else

            z=(median(x1-x2))/std(x1-x2);

        end



        r_val(i)=abs(z)/sqrt(minSubj);



        fprintf('%s vs %s: p=%.5f\n',...
            states{s1},states{s2},p);



    end



    %% ================= FDR =================


    p_fdr = mafdr(p_raw,...
        'BHFDR',true);



    %% ================= Final results =================


    fprintf('\n===== FINAL %s =====\n',...
        upper(band_name));



    for i=1:nPairs


        s1=pairs(i,1);

        s2=pairs(i,2);


        fprintf('%s vs %s | raw=%.5f | FDR=%.5f | r=%.3f\n',...
            states{s1},...
            states{s2},...
            p_raw(i),...
            p_fdr(i),...
            r_val(i));

    end



    %% ================= Plot =================


    figure('Color','w',...
           'Position',[200 200 750 550]);

    hold on;



    boxplot(X,...
        'Labels',stateLabels,...
        'Symbol','');



    boxes=findobj(gca,'Tag','Box');



    for i=1:length(boxes)


        patch(get(boxes(i),'XData'),...
              get(boxes(i),'YData'),...
              state_colors(length(boxes)-i+1,:),...
              'FaceAlpha',0.5,...
              'EdgeColor','k',...
              'LineWidth',1.2);

    end



    %% Individual trajectories


    for subj=1:minSubj


        plot(1:nStates,...
             X(subj,:),...
             '-o',...
             'Color',[0.7 0.7 0.7],...
             'LineWidth',1);

    end



    %% Mean


    plot(1:nStates,...
         mean(X,1),...
         'k-o',...
         'LineWidth',2,...
         'MarkerFaceColor','k');



    title([upper(band_name) ' band'],...
        'FontWeight','bold');



    ylabel('dwPLI');



    set(gca,...
        'FontSize',12,...
        'LineWidth',1.2);



    box off;



    %% ================= Significance =================


    ymax=max(X(:));

    ymin=min(X(:));

    yrange=ymax-ymin;



    h_gap=0.05*yrange;

    h_line=0.35*h_gap;



    y=ymax+h_gap;

    y_top=ymax;



    for i=1:nPairs


        if isnan(p_fdr(i)) || p_fdr(i)>=0.05

            continue;

        end



        x1=pairs(i,1);

        x2=pairs(i,2);



        if p_fdr(i)<0.001

            sig='***';

        elseif p_fdr(i)<0.01

            sig='**';

        else

            sig='*';

        end



        plot([x1 x1],...
             [y y+h_line],...
             'k','LineWidth',1.2);



        plot([x2 x2],...
             [y y+h_line],...
             'k','LineWidth',1.2);



        plot([x1 x2],...
             [y+h_line y+h_line],...
             'k','LineWidth',1.2);



        text(mean([x1 x2]),...
             y+1.5*h_gap,...
             sig,...
             'HorizontalAlignment','center',...
             'FontSize',12,...
             'FontWeight','bold');



        y_top=max(y_top,...
            y+1.5*h_gap);



        y=y+2*h_gap;



    end



    ylim([ymin,...
          y_top+0.8*h_gap]);



    %% ================= Save =================


    save_name=sprintf(...
        'dwPLI_FLC_HPC_%s.tiff',...
        band_name);



    print(gcf,...
        fullfile(data_dir,save_name),...
        '-dtiff',...
        '-r300');



    fprintf('Saved: %s\n',save_name);



end



disp('===== FLC-HPC dwPLI band statistics finished =====');