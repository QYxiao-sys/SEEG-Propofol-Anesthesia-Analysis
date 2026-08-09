%% ============================================================
% KLMI band-specific PAC extraction and statistics
%
% Description:
%   Extract frontal cortex (FLC) KLMI values in a predefined
%   frequency range, align subjects across states, perform
%   paired statistics with FDR correction, and generate boxplots.
%
% Frequency range:
%   Low frequency : 8-12 Hz
%   High frequency: 30-60 Hz
%
% Statistics:
%   Paired Wilcoxon signed-rank test
%   Benjamini-Hochberg FDR correction
%
% MATLAB version:
%   R2020b
%
% ============================================================


clc;
clear;
close all;



%% ========= Path =========

% KLMI result folder

dataDir = 'YOUR_KLMI_RESULT_FOLDER';


% Raw SEEG dataset folder

rawDir = 'YOUR_SEEG_DATA_ROOT';



stateNames  = {'MAIN','ANE','ETT','AWA','WARD'};

stateLabels = {'ANE','R-5','R+5','R+25','WARD'};



%% ========= Frequency band =========

lf_range = [8 12];

hf_range = [30 60];



%% ========= Colors =========

color.MAIN = [0.85 0.33 0.10];

color.ANE  = [0.93 0.69 0.13];

color.ETT  = [0.00 0.45 0.74];

color.AWA  = [0.47 0.67 0.19];

color.WARD = [0 0 0];


colors = [
    color.MAIN;
    color.ANE;
    color.ETT;
    color.AWA;
    color.WARD
];



%% ========= Storage =========

PAC_all  = cell(5,1);

Subj_all = cell(5,1);



%% =====================================================
% 1. Extract PAC + subject ID (FLC)
%% =====================================================


for s = 1:5


    state = stateNames{s};


    fprintf('\nProcessing %s...\n',state);



    file = fullfile(dataDir,...
        ['KLMI_CONTINUOUS_' state '.mat']);



    load(file);



    files = dir(fullfile(rawDir,state,'*.set'));

    names = {files.name};



    PAC_sub = [];

    subj_id = {};



    for k = 1:length(klmi_FLC)


        if isempty(klmi_FLC{k})

            continue;

        end



        img = klmi_FLC{k};



        lf_idx = find(...
            lf_lower >= lf_range(1) & ...
            lf_lower <= lf_range(2));



        hf_idx = find(...
            hf_lower >= hf_range(1) & ...
            hf_lower <= hf_range(2));



        val = mean(mean(img(lf_idx,hf_idx)));



        PAC_sub(end+1,1)=val;



        %% ===== Subject ID cleaning =====


        name = erase(names{k},'.set');


        name = regexprep(name,...
            '(MAIN|ANE|ETT|AWA|WARD)','');


        name = regexprep(name,...
            '[_\-]','');



        subj_id{end+1,1}=name;


    end



    PAC_all{s}=PAC_sub;

    Subj_all{s}=subj_id;


end



%% =====================================================
% 2. Subject alignment
%% =====================================================


common_subj = Subj_all{1};



for s = 2:5

    common_subj = intersect(...
        common_subj,...
        Subj_all{s});

end



fprintf('\n===== Common subjects: %d =====\n',...
    length(common_subj));



if isempty(common_subj)

    error('No common subjects found');

end



data_mat=zeros(length(common_subj),5);



for s=1:5


    for i=1:length(common_subj)


        idx=find(strcmp(...
            Subj_all{s},...
            common_subj{i}));


        data_mat(i,s)=PAC_all{s}(idx);


    end

end



%% =====================================================
% 3. Plot
%% =====================================================


figure('Color','w',...
    'Position',[200 200 700 550]);


hold on;



set(groot,...
    'defaultAxesFontName','Times New Roman');

set(groot,...
    'defaultTextFontName','Times New Roman');



boxplot(data_mat,...
    'Labels',stateLabels,...
    'Symbol','');



boxes=findobj(gca,'Tag','Box');


for j=1:length(boxes)


    patch(get(boxes(j),'XData'),...
          get(boxes(j),'YData'),...
          colors(length(boxes)-j+1,:),...
          'FaceAlpha',0.5,...
          'EdgeColor','k',...
          'LineWidth',1.2);


end



med=findobj(gca,'Tag','Median');


set(med,...
    'LineWidth',2,...
    'Color','k');



ymin=min(data_mat(:));

ymax=prctile(data_mat(:),95);


ylim([ymin ymax]);



ylabel('Phase-Amplitude Coupling Strength',...
    'FontWeight','bold');



title('Frontal cortex PAC (8-12 × 30-60 Hz)',...
    'FontSize',16,...
    'FontWeight','bold');



set(gca,...
    'FontSize',14,...
    'LineWidth',1.2);



ax=gca;

ax.YAxis.Exponent=-4;


box off



%% =====================================================
% 4. Statistics + significance annotation
%% =====================================================


nStates=size(data_mat,2);


pairs=nchoosek(1:nStates,2);


nComp=size(pairs,1);



pvals=zeros(nComp,1);



for i=1:nComp


    x=data_mat(:,pairs(i,1));

    y=data_mat(:,pairs(i,2));


    pvals(i)=signrank(x,y);


end



%% ===== Manual BH-FDR =====


[p_sorted,sort_idx]=sort(pvals);


V=length(pvals);


adj_p=zeros(size(pvals));


for i=1:V

    adj_p(i)=p_sorted(i)*V/i;

end



for i=V-1:-1:1

    adj_p(i)=min(adj_p(i),adj_p(i+1));

end



p_fdr=zeros(size(pvals));


p_fdr(sort_idx)=adj_p;



%% ===== Stars =====


yl=ylim;


y_base=yl(2);


y_step=(yl(2)-yl(1))*0.08;


count=0;



for i=1:nComp


    p=p_fdr(i);



    if p<0.001

        star='***';


    elseif p<0.01

        star='**';


    elseif p<0.05

        star='*';


    else

        continue;

    end



    count=count+1;



    x1=pairs(i,1);

    x2=pairs(i,2);



    y=y_base+count*y_step;



    plot([x1 x1 x2 x2],...
        [y y+y_step/4 y+y_step/4 y],...
        'k',...
        'LineWidth',1.5);



    text(mean([x1 x2]),...
        y+y_step*0.6,...
        star,...
        'HorizontalAlignment','center',...
        'FontSize',14,...
        'FontWeight','bold');


end



ylim([yl(1),...
      y_base+(count+1)*y_step]);



%% =====================================================
% 5. Save
%% =====================================================


outFig=fullfile(dataDir,...
    'FLC_PAC_boxplot.tif');


exportgraphics(gcf,...
    outFig,...
    'Resolution',300);



disp('===== DONE (FLC version) =====');