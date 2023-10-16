%% 20230727 GHR63-HC59 behav ana;
%% compute mean/SD of behavioral data
% GHRcolu = [5, 9, 11:14];% age,edu, TIVetc;
% HCcolu = [5, 9, 11:14];% age,sex,TIVetc;
% GHRcolu = [64, 65, 68, 72, 74, 76];% HVLT...;
% HCcolu = [36, 37, 40, 44, 46, 48];% HVLT...;
GHRcolu = [54, 55];
HCcolu = [26, 27];

GHRrow = [1:9, 11:63];% line 10 (GHR10181, abnormal LES)
HCrow = [1:11, 13:59];% line 12 (HC113, abnormal LES)
clear GHRmeansd;
clear HCmeansd;
for i = 1:length(GHRcolu)
    GHRmeansd(i,1) = nanmean(GHRdataMatchImageFDTIV1{GHRrow, GHRcolu(i)});
    GHRmeansd(i,2) = nanstd(GHRdataMatchImageFDTIV1{GHRrow, GHRcolu(i)});
    HCmeansd(i,1) = nanmean(HCdataMatchImageFDTIVspm1{HCrow, HCcolu(i)});
    HCmeansd(i,2) = nanstd(HCdataMatchImageFDTIVspm1{HCrow, HCcolu(i)});
end
% compute SIPS;
GHRcolu = [27, 33, 40, 45, 50, 51];% sips
for i = 1:length(GHRcolu)
    % cell2mat--cell->double
    GHRmeansd(i,1) = nanmean(GHRdataMatchImageFDTIV1{:, GHRcolu(i)});
    GHRmeansd(i,2) = nanstd(GHRdataMatchImageFDTIV1{:, GHRcolu(i)}); 
end

% t-test
% GHRcolu = [5, 9, 11:14];% age,edu, TIVetc;
% HCcolu = [5, 9, 11:14];% age,sex,TIVetc;
% GHRcolu = [64, 65, 68, 72, 74, 76];% HVLT...;
% HCcolu = [36, 37, 40, 44, 46, 48];% HVLT...;
GHRcolu = [54, 55];
HCcolu = [26, 27];
clear tt;
for j= 1: length(GHRcolu)
    [h,p,ci,stats] = ttest2(GHRdataMatchImageFDTIV1{GHRrow,GHRcolu(j)},...
        HCdataMatchImageFDTIVspm1{HCrow, HCcolu(j)});
    tt(j, 1) = stats.tstat; % ttest T value;
    tt(j, 2) = p; % ttest P value
end

%%anova WITH covariates
GHRcolu = [54, 55];
HCcolu = [26, 27];
Cov = [GHRdataMatchImageFDTIV1{GHRrow, [5, 4, 9]};...
    HCdataMatchImageFDTIVspm1{HCrow, [5, 4, 9]}];% age,sex,education;
clear fp;
for i = 1:length(GHRcolu)
    data2comp = [GHRdataMatchImageFDTIV1{GHRrow, GHRcolu(i)};...
        HCdataMatchImageFDTIVspm1{HCrow, HCcolu(i)}];
    grouplabel = [ones(length(GHRrow), 1);zeros(length(HCrow), 1)];%%%%%%%%%%%need revise;
    % age,sex,edu
    groups = [grouplabel, Cov];
    [p, tb1,stats]= anovan(data2comp, groups, 'continuous', [2, 4]);
    fp(i, 2) = p(1); % pvalue
    fp(i, 1) = tb1{2, 6}; % F value
end



%% %%%%%%%%%% subgroup GHR Behavioral analysis %%%%%%%%%%%;
%% subgroup the GHR by R_putamen volume-median of HC;
%% compute mean/SD of behavioral data
GHRcolu = [51, 27, 64, 65, 68, 72, 74,76];
GHRrowmax = GHRmaxindex78;
GHRrowmin = GHRminindex78;

% GHRcolu = [55];% LES,line 10 (GHR10181, abnormal LES)
% GHRrowmax = GHRmaxindex78;
% GHRrowmin = GHRminindex78([1:5, 7:37]);
clear GHRmeansd;
for i = 1:length(GHRcolu)
    % max subgroup;
    GHRmeansd(i,1) = nanmean(GHRdataMatchImageFDTIV1{GHRrowmax, GHRcolu(i)});
    GHRmeansd(i,2) = nanstd(GHRdataMatchImageFDTIV1{GHRrowmax, GHRcolu(i)});
    % min subgroup;
    GHRmeansd(i,3) = nanmean(GHRdataMatchImageFDTIV1{GHRrowmin, GHRcolu(i)});
    GHRmeansd(i,4) = nanstd(GHRdataMatchImageFDTIV1{GHRrowmin, GHRcolu(i)});
end

% compute hc mean/SD;
clear HCmeansd;
HCcolu = [36,37,40,44,46,48];
HCrow = [1:59];

% HCcolu = [26, 27];%line 12 (HC113, abnormal LES)
% HCrow = [1:11, 13:59];
for i = 1:length(HCcolu)
    HCmeansd(i,1) = nanmean(HCdataMatchImageFDTIVspm1{HCrow , HCcolu(i)});
    HCmeansd(i,2) = nanstd(HCdataMatchImageFDTIVspm1{HCrow , HCcolu(i)});  
end

%% anova WITH covariates--2023 between two subgroup;
colu1 =[51, 27]
row1 = GHRmaxindex78;% line 10 (GHR10181, abnormal LES)
row2 = GHRminindex78;% line 12 (HC113, abnormal LES)
% age, sex, education as COVs;
Cov = [GHRdataMatchImageFDTIV1{row1, [5, 4, 9]};...
    GHRdataMatchImageFDTIV1{row2, [5, 4, 9]}];
clear fp;
for i = 1:length(colu1)
    data2comp = [GHRdataMatchImageFDTIV1{row1, colu1(i)};...
        GHRdataMatchImageFDTIV1{row2, colu1(i)}];
    grouplabel = [ones(length(row1), 1);zeros(length(row2), 1)];%%%%%%%%%%%need revise;
    % age,sex,edu
    groups = [grouplabel, Cov];
    [p, tb1,stats]= anovan(data2comp, groups, 'continuous', [2, 4]);
    fp(i, 2) = p(1); % pvalue
    fp(i, 1) = tb1{2, 6}; % F value
end
%% anova WITH covariates--2023XIN(3 group + hoc-post);
% colu1 = [64, 65, 68, 72, 74,76];% HVLT, CPT, TMT. ..
% coluHC = [36, 37, 40, 44, 46, 48];
colu1 = [54, 55];% LES,line 10 (GHR10181, abnormal LES)
coluHC = [26, 27];%line 12 (HC113, abnormal LES)
row1 = GHRmaxindex78;
row2 = GHRminindex78([1:5, 7:37]);
rowHC = [1:11, 13:59];
Cov = [GHRdataMatchImageFDTIV1{row1, [5, 4, 9]};...
    GHRdataMatchImageFDTIV1{row2, [5, 4, 9]};...
    HCdataMatchImageFDTIVspm1{rowHC, [5, 4, 9]}];
clear fp tbl111;
for i = 1:length(colu1)
    data2comp = [GHRdataMatchImageFDTIV1{row1, colu1(i)};...
        GHRdataMatchImageFDTIV1{row2, colu1(i)};...
        HCdataMatchImageFDTIVspm1{rowHC, coluHC(i)}];
    grouplabel = [ones(length(row1), 1);2*ones(length(row2), 1);...
        3*ones(length(rowHC), 1)];%%%%%%%%%%%need revise;
    % age,sex,edu
    groups = [grouplabel, Cov];
    [p, tb1,stats]= anovan(data2comp, groups, 'continuous', [2, 4]);
    fp(i, 2) = p(1); % pvalue
    fp(i, 1) = tb1{2, 6}; % F value
    
    [results,~,~,gnames] = multcompare(stats,"Dimension",[1]);
    tbl111{i,1} = array2table(results,'VariableNames',...
    {'GroupA','GroupB' ,'LowerLimi','A_B','UpperLimit','Pvalue'});
end
