%% ============================================================
% Band-specific PSD statistical analysis
%
% Description:
%   Perform frequency-band level statistical analysis of ROI PSD
%   across five experimental states.
%
% Statistical methods:
%   - Friedman test
%   - Paired Wilcoxon signed-rank test
%   - Benjamini-Hochberg FDR correction
%   - Effect size r calculation
%
% Output:
%   TXT summary file
%   CSV statistical table
%
% MATLAB version:
%   R2020b
%
% ============================================================


clear;
clc;



disp('======================================')

disp('   PAIRED NONPARAMETRIC STATISTICS')

disp('======================================')



%% ===== ROI =====

ROI_list = {'PFC','HPC'};



%% ===== States =====

states      = {'MAIN','ANE','ETT','AWA','WARD'};

stateLabels = {'ANE','R-5','R+5','R+25','WARD'};



%% ===== Path =====

% Folder containing PSD_*.mat files

result_dir = 'YOUR_PSD_RESULT_FOLDER';



%% ===== TXT =====

txt_file = fullfile(result_dir,...
    'Band_Statistics_Paired_FDR.txt');


fid=fopen(txt_file,'w');



%% ===== CSV result container =====

results={};

result_count=0;



%% ===== Frequency bands =====

bands.delta=[1 4];

bands.theta=[4 8];

bands.alpha=[8 13];

bands.beta =[13 30];

bands.gamma=[30 80];


band_names=fieldnames(bands);



%% =====================================================
%% ================= MAIN LOOP ==========================
%% =====================================================


for r=1:length(ROI_list)


    roi=ROI_list{r};


    fprintf('\n================ ROI: %s ================\n',roi)


    fprintf(fid,...
        '\n================ ROI: %s ================\n',...
        roi);



    %% ===== load =====


    file=fullfile(result_dir,...
        ['PSD_' roi '.mat']);



    if ~exist(file,'file')


        warning('%s not found',file)

        continue


    end



    load(file);



    %% ===== check =====


    for k=1:length(states)


        if ~isfield(PSD,states{k})


            error('%s missing %s',...
                roi,states{k})


        end


    end



    %% ===== data =====


    data_MAIN=PSD.MAIN;

    data_ANE =PSD.ANE;

    data_ETT =PSD.ETT;

    data_AWA =PSD.AWA;

    data_WARD=PSD.WARD;



    %% =====================================================
    %% frequency loop
    %% =====================================================


    for b=1:length(band_names)



        band=band_names{b};

        range=bands.(band);



        fprintf('\n--- %s %.1f-%.1f Hz ---\n',...
            band,range(1),range(2))


        fprintf(fid,...
            '\n--- %s %.1f-%.1f Hz ---\n',...
            band,range(1),range(2));



        %% frequency index


        idx=freq>=range(1)&freq<=range(2);



        if sum(idx)==0

            continue

        end



        %% subject value


        MAIN=mean(data_MAIN(:,idx),2);

        ANE =mean(data_ANE(:,idx),2);

        ETT =mean(data_ETT(:,idx),2);

        AWA =mean(data_AWA(:,idx),2);

        WARD=mean(data_WARD(:,idx),2);



        band_data={MAIN,ANE,ETT,AWA,WARD};



        %% =================================================
        %% Friedman
        %% =================================================


        minN=min(cellfun(@length,band_data));


        data_mat=zeros(minN,length(states));


        for k=1:length(states)


            data_mat(:,k)=band_data{k}(1:minN);


        end



        p_friedman=friedman(data_mat,1,'off');



        fprintf(fid,...
            'Friedman p = %.5f\n',...
            p_friedman);



        %% =================================================
        %% Pairwise Wilcoxon
        %% =================================================


        nComp=nchoosek(length(states),2);



        p_vals=zeros(nComp,1);

        z_vals=zeros(nComp,1);

        r_vals=zeros(nComp,1);


        labels=cell(nComp,1);



        c=0;



        for i=1:length(states)-1


            for j=i+1:length(states)



                c=c+1;



                x=band_data{i};

                y=band_data{j};



                n=min(length(x),length(y));


                x=x(1:n);

                y=y(1:n);



                valid=~isnan(x)&~isnan(y);


                x=x(valid);

                y=y(valid);



                [p,~,stats_sr]=signrank(x,y);



                p_vals(c)=p;



                %% Z calculation


                d=x-y;


                d(d==0)=[];


                n_eff=length(d);



                if n_eff>0 && isfield(stats_sr,'signedrank')


                    W=stats_sr.signedrank;


                    mu=n_eff*(n_eff+1)/4;


                    sigma=sqrt(...
                        n_eff*(n_eff+1)*(2*n_eff+1)/24);



                    z_vals(c)=...
                        (W-mu)/sigma;



                    r_vals(c)=...
                        abs(z_vals(c))/sqrt(n_eff);



                else


                    z_vals(c)=NaN;

                    r_vals(c)=NaN;


                end



                labels{c}=...
                    [stateLabels{i} ' vs ' stateLabels{j}];


            end


        end



        %% ===== FDR =====


        p_fdr=mafdr(p_vals,...
            'BHFDR',true);



        %% =================================================
        %% save table
        %% =================================================


        for i=1:nComp



            result_count=result_count+1;



            results(result_count,:)={...
                roi,...
                band,...
                labels{i},...
                p_vals(i),...
                p_fdr(i),...
                r_vals(i)};



            if p_fdr(i)<0.001

                sig='***';

            elseif p_fdr(i)<0.01

                sig='**';

            elseif p_fdr(i)<0.05

                sig='*';

            else

                sig='ns';

            end



            fprintf('%s raw=%.4f FDR=%.4f r=%.3f %s\n',...
                labels{i},...
                p_vals(i),...
                p_fdr(i),...
                r_vals(i),...
                sig);



            fprintf(fid,...
                '%s raw=%.4f FDR=%.4f r=%.3f %s\n',...
                labels{i},...
                p_vals(i),...
                p_fdr(i),...
                r_vals(i),...
                sig);


        end


    end


end




%% =====================================================
%% Export CSV
%% =====================================================


result_table=cell2table(results,...
    'VariableNames',...
    {'ROI',...
    'frequency_band',...
    'Comparison',...
    'raw_p',...
    'FDR_p',...
    'Effect_size_r'});



csv_file=fullfile(result_dir,...
    'Band_Statistics_Paired_FDR_results.csv');



writetable(result_table,csv_file);



%% close


fclose(fid);



disp('======================================')

disp('FINISHED')

disp(['TXT : ' txt_file])

disp(['CSV : ' csv_file])

disp('======================================')