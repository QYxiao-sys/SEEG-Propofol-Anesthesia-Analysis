%% =====================================================
%%   FLC-HPC Directionality Index Statistics
%%
%%   Friedman test
%%   Pairwise Wilcoxon signed-rank test
%%   FDR correction
%%   Effect size r
%%
%%   Export TXT + CSV
%%
%%   Required data structure:
%%
%%   data/
%%      MAIN/
%%      ANE/
%%      ETT/
%%      AWA/
%%      WARD/
%%
%%   Each folder contains:
%%      MVGC_*.mat
%%
%%   Each MAT file contains:
%%      result.GC
%%      result.freq
%%
%% =====================================================


clear; clc; close all;



%% =====================================================
%% Relative paths
%% =====================================================


codeDir = fileparts(mfilename('fullpath'));

dataDir = fullfile(codeDir,'..','data');

resultDir = fullfile(codeDir,'..','results');


if ~exist(resultDir,'dir')
    mkdir(resultDir);
end



%% ================= ROI =================


ROI = 'FLC-HPC';



%% ================= States =================


states = {'MAIN','ANE','ETT','AWA','WARD'};

stateLabels = {'ANE','R-10','R+5','R+25','WARD'};

nStates = length(states);



%% ================= TXT output =================


txt_file = fullfile(resultDir,...
    'FLC_HPC_DI_statistics_FDR.txt');


fid = fopen(txt_file,'w');



%% ================= CSV container =================


results = {};

result_count = 0;



%% ================= Frequency bands =================


bands = { ...
    'delta',[1 4]; ...
    'theta',[4 8]; ...
    'alpha',[8 13]; ...
    'beta',[13 30]; ...
    'gamma',[30 100]};



%% =====================================================
%% Main analysis
%% =====================================================


for b = 1:size(bands,1)


    band_name = bands{b,1};

    band_range = bands{b,2};



    fprintf('\n========== %s ==========\n',...
        upper(band_name));


    fprintf(fid,...
        '\n========== %s ==========\n',...
        upper(band_name));



    %% =========================
    %% Extract DI for each state
    %% =========================


    DI_data = cell(nStates,1);



    for s = 1:nStates


        stateDir = fullfile(dataDir,states{s});


        files = dir(fullfile(stateDir,...
            'MVGC_*.mat'));



        allDI = [];



        for fIdx = 1:length(files)



            try


                load(fullfile(stateDir,...
                    files(fIdx).name),...
                    'result');



                GC = result.GC;

                freq = result.freq;



                %% Frequency selection

                freq_idx = freq>=1 & freq<=100;


                freq_plot = freq(freq_idx);



                %% =====================================
                % Channel order:
                %
                % 1 = FLC
                % 2 = HPC
                %
                % GC(1,2):
                % FLC → HPC
                %
                % GC(2,1):
                % HPC → FLC
                %% =====================================


                gc_FLC2HPC = squeeze(GC(1,2,:));

                gc_HPC2FLC = squeeze(GC(2,1,:));



                gc_FLC2HPC = gc_FLC2HPC(freq_idx);

                gc_HPC2FLC = gc_HPC2FLC(freq_idx);



                %% =========================
                % Directionality Index
                %
                % Positive:
                % HPC → FLC dominance
                %
                % Negative:
                % FLC → HPC dominance
                %% =========================


                DI = (gc_HPC2FLC - gc_FLC2HPC) ./ ...
                     (gc_HPC2FLC + gc_FLC2HPC + eps);



                %% Frequency band extraction


                band_idx = freq_plot>=band_range(1) & ...
                           freq_plot<=band_range(2);



                DI_band = mean(DI(band_idx));



                allDI = [allDI; DI_band];



            catch

                warning('Failed: %s',...
                    files(fIdx).name);

                continue

            end


        end



        DI_data{s}=allDI;



    end



    %% =========================
    %% Arrange paired data
    %% =========================


    minSubj = min(cellfun(@length,DI_data));


    if isempty(minSubj) || minSubj==0

        warning('No valid subjects in %s',band_name);

        continue

    end



    X=zeros(minSubj,nStates);



    for s=1:nStates

        X(:,s)=DI_data{s}(1:minSubj);

    end




    %% =========================
    %% Friedman test
    %% =========================


    p_friedman = friedman(X,1,'off');


    fprintf('Friedman p = %.6f\n',...
        p_friedman);


    fprintf(fid,...
        'Friedman p = %.6f\n',...
        p_friedman);



    %% =========================
    %% Pairwise Wilcoxon
    %% =========================


    pairs = nchoosek(1:nStates,2);

    nPairs = size(pairs,1);


    p_raw=zeros(nPairs,1);

    r_val=zeros(nPairs,1);

    comparisons=cell(nPairs,1);



    for i=1:nPairs



        s1=pairs(i,1);

        s2=pairs(i,2);



        x1=X(:,s1);

        x2=X(:,s2);



        [p,~,stats]=signrank(x1,x2);


        p_raw(i)=p;



        %% Effect size r

        d=x1-x2;

        d(d==0)=[];

        n_eff=length(d);



        if n_eff>0 && isfield(stats,'signedrank')


            W=stats.signedrank;


            mu=n_eff*(n_eff+1)/4;


            sigma=sqrt(...
                n_eff*(n_eff+1)*(2*n_eff+1)/24);



            Z=(W-mu)/sigma;


            r_val(i)=abs(Z)/sqrt(n_eff);



        else

            r_val(i)=NaN;

        end



        comparisons{i}=...
            [stateLabels{s1} ' vs ' stateLabels{s2}];



    end




    %% =========================
    %% FDR correction
    %% =========================


    p_fdr = mafdr(p_raw,...
        'BHFDR',true);




    %% =========================
    %% Save results
    %% =========================


    for i=1:nPairs



        result_count=result_count+1;



        results(result_count,:)={...
            ROI,...
            band_name,...
            comparisons{i},...
            p_raw(i),...
            p_fdr(i),...
            r_val(i)};



        if p_fdr(i)<0.001

            sig='***';

        elseif p_fdr(i)<0.01

            sig='**';

        elseif p_fdr(i)<0.05

            sig='*';

        else

            sig='ns';

        end



        fprintf('%s | raw %.5f | FDR %.5f | r %.3f %s\n',...
            comparisons{i},...
            p_raw(i),...
            p_fdr(i),...
            r_val(i),...
            sig);



        fprintf(fid,...
            '%s | raw %.5f | FDR %.5f | r %.3f %s\n',...
            comparisons{i},...
            p_raw(i),...
            p_fdr(i),...
            r_val(i),...
            sig);


    end


end




%% =====================================================
%% Export CSV
%% =====================================================


result_table = cell2table(results,...
    'VariableNames',...
    {'ROI',...
     'frequency_band',...
     'Comparison',...
     'raw_p',...
     'FDR_p',...
     'Effect_size_r'});



csv_file = fullfile(resultDir,...
    'FLC_HPC_DI_statistics_FDR_results.csv');



writetable(result_table,csv_file);



fclose(fid);



disp('======================================')
disp('FLC-HPC DI statistics completed')
disp(['TXT: ',txt_file])
disp(['CSV: ',csv_file])
disp('======================================')