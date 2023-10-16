%% 1. select the T1 image for GHR63 HC59, add computed TIV data to the Behavioral data;% UFR/HC
dataA = importdata('TIV.txt');
dataAta = array2table(dataA(:, [1:4]),'VariableNames', {'TIV', 'GM',...
    'WM', 'CSF'});
HCdataMatchImageFDTIVspm1 = [HCdataMatchImageFDTIVspm(:, [1:10]), dataAta,...
    HCdataMatchImageFDTIVspm(:, [15: end])];

%% compute the volume of subarea of striatum
% step1: obtain the ROI seed(subarea of striatum in AAL3);
path = ['~/3ndComputVBMcat12modulated'];
mask1 = zeros(113,137,113);
V1 = spm_vol([path, '/slice_AAL3v1.nii']);% 1.5*1.5*1.5;
T1 = spm_read_vols(V1);
for i = 157:158
    maskX = mask1+ (T1==i);
    sliceN = ['slice_AAL3striatum', num2str(i), 'mask.nii']
    V1.fname =  [path, filesep, sliceN];
    mask2 = spm_write_vol(V1, maskX);
end
% step2:extract volume of the subarea of striatum 
workpath = ['~/3ndComputVBMcat12modulated'];
groupP = [workpath, '/GHR63T1/mri'];%%%%%%%GHR63/HC59;
subfolder = dir([groupP, '/s6mwp1sub-*_T1w.nii']);

volume = zeros(113,137,113); 
roiVol = zeros(length(subfolder), 7);

modicell = {'slice_AAL3striatum75mask.nii', 'slice_AAL3striatum76mask.nii',...
        'slice_AAL3striatum77mask.nii', 'slice_AAL3striatum78mask.nii',...
        'slice_AAL3striatum157mask.nii', 'slice_AAL3striatum158mask.nii'};
clear roiVol;
for i = 1:length(subfolder)
    splitcell = split(subfolder(i).name, '_');
    split1cell = split(splitcell{1}, '-');
    roiVol(i, 1) = str2num(split1cell{2});
    Vsub = spm_vol([groupP, '/', subfolder(i).name]);
    volume1 = spm_read_vols(Vsub);
    
    for roii=1:length(modicell)
        clear plus;
        Vroi = spm_vol([workpath, '/', modicell{roii}]);
        volume = spm_read_vols(Vroi);
        plus = volume1.*volume;
        total = sum(plus(:));
        roiVol(i, roii+1) = total* 1.5^3/10^3;% mm3 --> ml
    end
end
roiVoltable = array2table(roiVol, 'VariableNames',...
    {'ID', 'cau75volml', 'cau76volml', 'put77volml',...
    'put78colml', 'NAcc157volml', 'NAcc158volml'});

%% mean and SD of striatum volumes(or GM WM CSF) for GHR63-HC59;
colu = [2: 7];
clear meansd;
for i = 1:length(colu)
    % group;
    meansd(i,1) = nanmean(GHR63subareaStriatumVOL1{:, colu(i)});
    meansd(i,2) = nanstd(GHR63subareaStriatumVOL1{:, colu(i)});
    meansd(i,3) = nanmean(HC59subareaStriatumVOL1{:, colu(i)});%%%%HC
    meansd(i,4) = nanstd(HC59subareaStriatumVOL1{:, colu(i)});
end
%% anova WITH covariates--2023 between GHR and HC volume;
colu1 =[2: 7];
row1 = [1:63];
row2 = [1:59];
Cov = [GHRdataMatchImageFDTIV1{row1, [5, 4, 9, 11]};...
    HCdataMatchImageFDTIVspm1{row2, [5, 4, 9,11]}];
clear fp;
for i = 1:length(colu1)
    data2comp = [GHR63subareaStriatumVOL1{row1, colu1(i)};...
       HC59subareaStriatumVOL1{row2, colu1(i)}];
    grouplabel = [ones(length(row1), 1);zeros(length(row2), 1)];%%%%%%%%%%%need revise;
    % age,sex,edu, TIV;
    groups = [grouplabel, Cov];
    [p, tb1,stats]= anovan(data2comp, groups, 'continuous', [2, 4, 5]);
    fp(i, 2) = p(1); % pvalue
    fp(i, 1) = tb1{2, 6}; % F value
end
%% correaltion between striatum subarea with SIPS-GHR63;
rowINDEX  = {[1:63]};% 1*4 cell;
[r, c] = size(rowINDEX);

col = [27, 51];% col index of SIPS;
% col = [55]; % LESN
% col= [77:85]; % congitive
striROIco = [2:7];% column index of striatum subarea vol;
clear RPGHR;
% for each group;
for groupi = 1: c 
    clear Cov;
    % Age,sex,edu, TIV;
    Cov = GHRdataMatchImageFDTIV1{rowINDEX{1, groupi}, [5, 4, 9, 11]};
    
    for  i =1: length(striROIco) % for each subarea of striatum
        for j = 1:length(col)
            colu = j*2;
            [R, P] = partialcorr(GHR63subareaStriatumVOL1{rowINDEX{1, groupi}, striROIco(i)},...
                GHRdataMatchImageFDTIV1{rowINDEX{1, groupi}, col(j)},...
                Cov, 'Type', 'Spearman', 'Rows', 'Complete');
            RPGHR(6*(groupi-1)+i, colu-1) = R; % need revise;
            RPGHR(6*(groupi-1)+i, colu) = P; % need revise
        end
    end
end
fdr = mafdr(RPGHR([1: 6], 2), 'BHFDR', true);
%% %% sctter the striatum subarea with SIPS-GHR63; (residual regressed age,sex,edu,TIV);
% volume-SIPS;
indexGHR = [1:63];
row = [2:3];%Xdata
variab = {'Lcau', 'Rcau'};
col = [27];
variabr = {'SIPStotal'};

figure('color', 'w')
for ri = 1:length(row)
    for c = 1:length(col)
        subplot(length(row), length(col), (ri-1)+c);
        
        Xdata = GHR63subareaStriatumVOL1{indexGHR, row(ri)}; %%%%VOL, need revise;
        Ydata = GHRdataMatchImageFDTIV1{indexGHR, col(c)}; %%%% SIPS/GAF, need revise;
        
        X = [ones(length(indexGHR), 1) GHRdataMatchImageFDTIV1{indexGHR, [5, 4, 9, 11]}]; %%%%% age,sex,edu
        [b,bint,r] = regress(Xdata, X); %
        [b1,bint1,r1] = regress(Ydata, X);%
        Xres(:, ri) = r;
        Yres(:, c) = r1;
        
        theta = glmfit(r, r1);
        yCalc1 = theta(2)*r +theta(1);
        plot(r(:, 1), r1(:, 1), 'Color', sixteen2ten('#B6825F')/255,...
            'Marker','o', 'Markersize',6, 'LineStyle', 'none', 'LineWidth', 1);
        hold on;
        
        plot(r(:, 1), yCalc1, 'Color', sixteen2ten('#FACB66')/255, 'LineWidth', 2);
        
        ylabel(['Res ', variabr{c}], 'FontSize', 8, 'FontWeight', 'bold');
        xlabel(['Res ', variab{ri}], 'FontSize', 8, 'FontWeight', 'bold');
        hold on
    end
end
% Original volume-SIPS;
indexGHR = [1:63];
row = [2:3];%Xdata
variab = {'Lcau', 'Rcau'};
col = [27];
variabr = {'SIPStotal'};

figure('color', 'w')
for ri = 1:length(row)
    for c = 1:length(col)
        subplot(length(row), length(col), (ri-1)+c);
        
        Xdata = GHR63subareaStriatumVOL1{indexGHR, row(ri)}; %%%%VOL, need revise;
        Ydata = GHRdataMatchImageFDTIV1{indexGHR, col(c)}; %%%% SIPS/GAF, need revise;
        
        theta = glmfit(Xdata, Ydata);
        yCalc1 = theta(2)*Xdata +theta(1);
        plot(Xdata(:, 1), Ydata(:, 1), 'Color', sixteen2ten('#B6825F')/255,...
            'Marker','*', 'Markersize',6, 'LineStyle', 'none', 'LineWidth', 1);
        hold on;
        
        plot(Xdata(:, 1), yCalc1, 'Color', sixteen2ten('#FACB66')/255, 'LineWidth', 2);
        
        ylabel(['Original ', variabr{c}], 'FontSize', 8, 'FontWeight', 'bold');
        xlabel(['Original ', variab{ri}], 'FontSize', 8, 'FontWeight', 'bold');
        hold on
    end
end

 %% 202307: GHR63: correlation between striatum subarea with LESN;
rowINDEX  = {[1:9, 11:63]};%lin10 of GHR63, abnormal LES;
[r, c] = size(rowINDEX);
col = [55];
striROIco = [2:7];% column index of striatum subarea vol;
clear RPGHR;
% for each group;
for groupi = 1: c 
    clear Cov;
    % Age,sex,edu, TIV;
    Cov = GHRdataMatchImageFDTIV1{rowINDEX{1, groupi}, [5, 4, 9, 11]};
    
    for  i =1: length(striROIco) % for each subarea of striatum
        for j = 1:length(col)
            colu = j*2;
            [R, P] = partialcorr(GHR63subareaStriatumVOL1{rowINDEX{1, groupi}, striROIco(i)},...
                GHRdataMatchImageFDTIV1{rowINDEX{1, groupi}, col(j)},...
                Cov, 'Type', 'Spearman', 'Rows', 'Complete');
            RPGHR(8*(groupi-1)+i, colu-1) = R; % need revise;
            RPGHR(8*(groupi-1)+i, colu) = P; % need revise
        end
    end
end

%% %%%%subgroup the GHR by striatum(including putamen...) volume-median from HC;
%% compute the mean(SD) striatum(including putamen...) volume of HC59;
colu = 7;
meansd(1,1) = nanmean(HC59subareaStriatumVOL1{:, colu});
meansd(1,2) = nanstd(HC59subareaStriatumVOL1{:, colu});
meansd(1,3) = median(HC59subareaStriatumVOL1{:, colu});
% GHR63
meansd(2,1) = nanmean(GHR63subareaStriatumVOL1{:, colu});
meansd(2,2) = nanstd(GHR63subareaStriatumVOL1{:, colu});
meansd(2,3) = median(GHR63subareaStriatumVOL1{:, colu});
%% %%%%%%%%%%%%%%%1. separate the GHR by cutoff defined as the median of volume in HC;
clear num;
cutof =0.6681; % vol median;
col = 7; % col=2 ---L_caudate;

num(1, 1) = sum(GHR63subareaStriatumVOL1{:, col}>cutof)
num(2, 1) = sum(GHR63subareaStriatumVOL1{:, col}<=cutof)

num(1, 2) = sum(HC59subareaStriatumVOL1{:, col}>cutof)
num(2, 2) = sum(HC59subareaStriatumVOL1{:, col}<=cutof)

%% 1) subgroup the GHR by caudate/putamen/ NAcc  volume-median of HCs;
% load GHR and HC volume table;
load('GHR63subareaStriatumVOL1.mat');
load('HC59subareaStriatumVOL1.mat');
% (1) putamen-median as cutoff, separately;
cutof = 0.6681;% vol median mean;
col = 7;% col=2 ---L_caudate;

GHRmaxindex158 = find(GHR63subareaStriatumVOL1{:, col}>cutof);
save GHRmaxHCmedian_index158.mat GHRmaxindex158; 
GHRminindex158 = find(GHR63subareaStriatumVOL1{:, col}<=cutof);
save GHRminHCmedian_index158.mat GHRminindex158; 

HCmaxindex158 = find(HC59subareaStriatumVOL1{:, col}>cutof);
save HCmaxHCmedian_index158.mat HCmaxindex158; 
HCminindex158 = find(HC59subareaStriatumVOL1{:, col}<=cutof);
save HCminHCmedian_index158.mat HCminindex158; 

%% mean(SD) of subarea subgroup
colu = [5, 9, 11:14];
row1 = GHRmaxindex78;
row2 = GHRminindex78;
rowHC = [1:59];
clear meansd;
for i = 1:length(colu)
    % group;
    meansd(i,1) = nanmean(GHRdataMatchImageFDTIV1{row1, colu(i)});
    meansd(i,2) = nanstd(GHRdataMatchImageFDTIV1{row1, colu(i)});
    meansd(i,3) = nanmean(GHRdataMatchImageFDTIV1{row2, colu(i)});
    meansd(i,4) = nanstd(GHRdataMatchImageFDTIV1{row2, colu(i)});
    meansd(i,5) = nanmean(HCdataMatchImageFDTIVspm1{rowHC, colu(i)});%%%%HC
    meansd(i,6) = nanstd(HCdataMatchImageFDTIVspm1{rowHC, colu(i)});
end
%% anova WITH covariates--2023(3 group + hoc-post);
%  behavioral variable, symptom, cognitive,,,
rowHC = [1:59];
colu1 =[12:14];
row1 = GHRmaxindex78;%
row2 = GHRminindex78;
Cov = [GHRdataMatchImageFDTIV1{row1, [5, 4, 9, 11]};...
    GHRdataMatchImageFDTIV1{row2, [5, 4, 9, 11]};...
    HCdataMatchImageFDTIVspm1{rowHC, [5, 4, 9, 11]}];
clear fp;
clear tbl11;
for i = 1:length(colu1)
    data2comp = [GHRdataMatchImageFDTIV1{row1, colu1(i)};...
        GHRdataMatchImageFDTIV1{row2, colu1(i)};...
        HCdataMatchImageFDTIVspm1{rowHC, colu1(i)}];
    grouplabel = [ones(length(row1), 1);2*ones(length(row2), 1);...
        3*ones(length(rowHC), 1)];%%%%%%%%%%%need revise;
    % age,sex,edu,tiv
    groups = [grouplabel, Cov];
    [p, tb1,stats]= anovan(data2comp, groups, 'continuous', [2, 4, 5]);
    fp(i, 2) = p(1); % pvalue
    fp(i, 1) = tb1{2, 6}; % F value
    
    [results,~,~,gnames] = multcompare(stats,"Dimension",[1]);
    tbl11{i, 1} = array2table(results,'VariableNames',...
    {'GroupA','GroupB' ,'LowerLimi','A_B','UpperLimit','Pvalue'});
end

%% (subgroup of striatum volume) anova WITH covariates--2023(3 group + hoc-post);
rowHC = [1:59];
colu1 =[2:7];
row1 = GHRmaxindex78;%
row2 = GHRminindex78;
Cov = [GHRdataMatchImageFDTIV1{row1, [5, 4, 9, 11]};...
    GHRdataMatchImageFDTIV1{row2, [5, 4, 9, 11]};...
    HCdataMatchImageFDTIVspm1{rowHC, [5, 4, 9, 11]}];
clear fp;
clear tbl11;
for i = 1:length(colu1)
    data2comp = [GHR63subareaStriatumVOL1{row1, colu1(i)};...
        GHR63subareaStriatumVOL1{row2, colu1(i)};...
        HC59subareaStriatumVOL1{rowHC, colu1(i)}];
    grouplabel = [ones(length(row1), 1);2*ones(length(row2), 1);...
        3*ones(length(rowHC), 1)];%%%%%%%%%%%need revise;
    % age,sex,edu,tiv
    groups = [grouplabel, Cov];
    [p, tb1,stats]= anovan(data2comp, groups, 'continuous', [2, 4, 5]);
    fp(i, 2) = p(1); % pvalue
    fp(i, 1) = tb1{2, 6}; % F value
    
    [results,~,~,gnames] = multcompare(stats,"Dimension",[1]);
    tbl11{i, 1} = array2table(results,'VariableNames',...
    {'GroupA','GroupB' ,'LowerLimi','A_B','UpperLimit','Pvalue'});
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) (subgroup1-HCs: comparing striatum volume) anova WITH covariates--20230730
clear Covar GHRrow datavol1 datavol2 meansdT1;
GHRrow = GHRmaxindex158;
HCrow = [1: 59];
GHRcolu = [2: 7];
HCcolu = [2: 7];
Covar = [GHRdataMatchImageFDTIV1{GHRrow, [5, 4, 9, 11]};HCdataMatchImageFDTIVspm1{HCrow, [5, 4, 9, 11]}];
clear fp1;
for i = 1:length(GHRcolu)
    data2comp = [GHR63subareaStriatumVOL1{GHRrow, GHRcolu(i)};...
        HC59subareaStriatumVOL1{HCrow, HCcolu(i)}];
    grouplabel = [ones(length(GHRrow), 1);zeros(length(HCrow), 1)];%%%%%%%%%%%need revise;
    % age,sex,edu
    groups = [grouplabel, Covar];
    [p, tb1,stats]= anovan(data2comp, groups, 'continuous', [2, 4, 5]);
    fp1(i, 2) = p(1); % pvalue
    fp1(i, 1) = tb1{2, 6}; % F value
    
    % mean,sd, ttest,p
    datavol1 = GHR63subareaStriatumVOL1{GHRrow, GHRcolu(i)};% GHR subgroup;
    datavol2= HC59subareaStriatumVOL1{HCrow, HCcolu(i)};% HC subgroup;
    
    meansdT1(i,1) = nanmean(datavol1);
    meansdT1(i,2) = nanstd(datavol1);
    meansdT1(i,3) = nanmean(datavol2);
    meansdT1(i,4) = nanstd(datavol2);
    
      [h,p,ci,stats] = ttest2(datavol1, datavol2);
    meansdT1(i, 5) = stats.tstat; % ttest T value;
    meansdT1(i, 6) = p; % ttest P value
end

%%(subgroup1-HCs: comparing striatum volume)
clear Covar GHRrow datavol1 datavol2 meansdT2;
GHRrow = GHRminindex158;
HCrow = [1: 59];
GHRcolu = [2: 7];
HCcolu = [2: 7];
Covar = [GHRdataMatchImageFDTIV1{GHRrow, [5, 4, 9, 11]}; HCdataMatchImageFDTIVspm1{HCrow, [5, 4, 9, 11]}];
clear fp2;
for i = 1:length(GHRcolu)
    data2comp = [GHR63subareaStriatumVOL1{GHRrow, GHRcolu(i)};...
        HC59subareaStriatumVOL1{HCrow, HCcolu(i)}];
    grouplabel = [ones(length(GHRrow), 1);zeros(length(HCrow), 1)];%%%%%%%%%%%need revise;
    % age,sex,edu
    groups = [grouplabel, Covar];
    [p, tb1,stats]= anovan(data2comp, groups, 'continuous', [2, 4, 5]);
    fp2(i, 2) = p(1); % pvalue
    fp2(i, 1) = tb1{2, 6}; % F value
    
    % mean,sd, ttest,p
    datavol1 = GHR63subareaStriatumVOL1{GHRrow, GHRcolu(i)};% GHR subgroup;
    datavol2= HC59subareaStriatumVOL1{HCrow, HCcolu(i)};% HC subgroup;
    
    meansdT2(i,1) = nanmean(datavol1);
    meansdT2(i,2) = nanstd(datavol1);
    meansdT2(i,3) = nanmean(datavol2);
    meansdT2(i,4) = nanstd(datavol2);
    
      [h,p,ci,stats] = ttest2(datavol1, datavol2);
    meansdT2(i, 5) = stats.tstat; % ttest T value;
    meansdT2(i, 6) = p; % ttest P value
end

%% 3)correlation of subgroup striatum(subarea) vol and SIPS,,,LES;
 rowINDEX  = {GHRmaxindex78([1:20, 22:26]), GHRminindex78([1:5, 7:37])};% 1*4 cell;(GHR10181, abnormal LES)
[r, c] = size(rowINDEX);

col = [55];%LES
% col = [27, 51]; % sips, GAF;

striROIco = [2: 7];% column index of striatum subarea vol;
N = length(striROIco);

clear RPGHR;
% for each group;
for groupi = 1: c 
    clear Cov;
    % Age,sex,edu, TIV;
    Cov = GHRdataMatchImageFDTIV1{rowINDEX{1, groupi}, [5, 4, 9, 11]};
    
    for  i =1: length(striROIco) % for each subarea of striatum
        for j = 1:length(col)
            colu = j*2;
            [R, P] = partialcorr(GHR63subareaStriatumVOL1{rowINDEX{1, groupi}, striROIco(i)},...
                GHRdataMatchImageFDTIV1{rowINDEX{1, groupi}, col(j)},...
                Cov, 'Type', 'Spearman', 'Rows', 'Complete');
            RPGHR(N*(groupi-1)+i, colu-1) = R; 
            RPGHR(N*(groupi-1)+i, colu) = P; 
        end
    end
end

%% other separation: correlation of subgroup striatum(subarea) vol and SIPS,,,LES;
rowINDEX  = {GHRmaxindex158, GHRminindex158([1:6, 8:37])};
[r, c] = size(rowINDEX);
 col = [55];%LES
striROIco = [2:7];% column index of striatum subarea vol;
N = length(striROIco);

clear RPGHR;
% for each group;
for groupi = 1: c 
    clear Cov;
    % Age,sex,edu, TIV;
    Cov = GHRdataMatchImageFDTIV1{rowINDEX{1, groupi}, [5, 4, 9, 11]};
    
    for  i =1: length(striROIco) % for each subarea of striatum
        for j = 1:length(col)
            colu = j*2;
            [R, P] = partialcorr(GHR63subareaStriatumVOL1{rowINDEX{1, groupi}, striROIco(i)},...
                GHRdataMatchImageFDTIV1{rowINDEX{1, groupi}, col(j)},...
                Cov, 'Type', 'Spearman', 'Rows', 'Complete');
            RPGHR(N*(groupi-1)+i, colu-1) = R; 
            RPGHR(N*(groupi-1)+i, colu) = P; 
        end
    end
end
%% scatter sign two subgroups (residual regressed age,sex,edu,TIV);
% volume-SIPS;
indexGHRmin = GHRminindex78([1:5, 7:37]);

col = [2, 3, 4, 6];%Xdata
variab = {'Lcau', 'Rcauda', 'Lputa','LNacc'};
row = [55];% Ydata
variabr = {'LESN'};

figure('color', 'w')
for ri = 1:length(row)
    for c = 1:length(col)
        subplot(2, 2, c);      
        % min group;
        Xdata1 = GHR63subareaStriatumVOL1{indexGHRmin, col(c)}; %%%%VOL, need revise;
        Ydata1 = GHRdataMatchImageFDTIV1{indexGHRmin, row(ri)}; %%%% SIPS/GAF, need revise;
        
        X1 = [ones(length(indexGHRmin), 1) GHRdataMatchImageFDTIV1{indexGHRmin, [5, 4, 9, 11]}]; %%%%% age,sex,edu
        [b,bint,rx] = regress(Xdata1, X1); %
        [b1,bint1,r1x] = regress(Ydata1, X1);%
        
        thetax = glmfit(rx, r1x);
        yCalc11 = thetax(2)*rx +thetax(1);
        plot(rx(:, 1), r1x(:, 1), 'Color', sixteen2ten('#718BC1')/255,...
            'Marker','o', 'Markersize',8, 'LineStyle', 'none', 'LineWidth', 1);
        hold on;
        
        plot(rx(:, 1), yCalc11, 'Color', sixteen2ten('#586B9C')/255, 'LineWidth', 2); 
        ylabel(['Res ', variabr{ri}], 'FontSize', 8, 'FontWeight', 'bold');
        xlabel(['Res ', variab{c}], 'FontSize', 8, 'FontWeight', 'bold');
    end
end

%% scatter two subgroups (original score without cov);
%  GHR78 Max: vol ---SIPS/LESN;
indexGHRmin = GHRmaxindex78;

col = [5];%Xdata
variab = {'Rputa'};
row = [27];% Ydata
variabr = {'SIPStot'};

figure('color', 'w')
for r = 1: length(row)
    for c = 1:length(col)
        %subplot(2, 2, c);
        % min group;
        Xdata1 = GHR63subareaStriatumVOL1{indexGHRmin, col(c)}; %%%%VOL, need revise;
        Ydata1 = GHRdataMatchImageFDTIV1{indexGHRmin, row(r)}; %%%% SIPS/GAF, need revise;
        
        thetax = glmfit(Xdata1, Ydata1);
        yCalc11 = thetax(2)*Xdata1 +thetax(1);
        plot(Xdata1(:, 1), Ydata1(:, 1), 'Color', sixteen2ten('#718BC1')/255,...
            'Marker','*', 'Markersize',6, 'LineStyle', 'none', 'LineWidth', 1);
        hold on;
        
        plot(Xdata1(:, 1), yCalc11, 'Color', sixteen2ten('#586B9C')/255, 'LineWidth', 2);
        ylabel(['Orig ', variabr{ri}], 'FontSize', 8, 'FontWeight', 'bold');
        xlabel(['Orig ', variab{c}], 'FontSize', 8, 'FontWeight', 'bold');
    end
end
%% scatter sign two subgroups (residual regressed age,sex,edu,TIV);
% volume-SIPS;
indexGHRmin = GHRmaxindex78;

col = [5];%Xdata
variab = {'Rputa'};
row = [27];% Ydata
variabr = {'SIPStot'};

figure('color', 'w')
for ri = 1:length(row)
    for c = 1:length(col)
        % min group;
        Xdata1 = GHR63subareaStriatumVOL1{indexGHRmin, col(c)}; %%%%VOL, need revise;
        Ydata1 = GHRdataMatchImageFDTIV1{indexGHRmin, row(ri)}; %%%% SIPS/GAF, need revise;
        
        X1 = [ones(length(indexGHRmin), 1) GHRdataMatchImageFDTIV1{indexGHRmin, [5, 4, 9, 11]}]; %%%%% age,sex,edu
        [b,bint,rx] = regress(Xdata1, X1); %
        [b1,bint1,r1x] = regress(Ydata1, X1);%
        
        thetax = glmfit(rx, r1x);
        yCalc11 = thetax(2)*rx +thetax(1);
        plot(rx(:, 1), r1x(:, 1), 'Color', sixteen2ten('#718BC1')/255,...
            'Marker','o', 'Markersize',8, 'LineStyle', 'none', 'LineWidth', 1);
        hold on;
        
        plot(rx(:, 1), yCalc11, 'Color', sixteen2ten('#586B9C')/255, 'LineWidth', 2); 
        ylabel(['Res ', variabr{ri}], 'FontSize', 8, 'FontWeight', 'bold');
        xlabel(['Res ', variab{c}], 'FontSize', 8, 'FontWeight', 'bold');
    end
end

%% %%%%%% GHR63-HC59 VBM Analysis;%%%%%%%%%%%%%%%%%%%%%
%% %%% volume analysis between each GHR63 group and HC59 %%%
%% 1. to extract the cluster volume in the group of GHR63 and HC59, separately;
%%%%%%%%%%%%%%% GHR81 HC59 VBM 
WKD = ['/public/home/Fengnana/DATA210725/SCZ2015XYMRS/PROdcmniixBIDSnoMRS/BIDSNEWfmriprep/3ndComputVBMcat12modulated/ghr63HC59volume2ndanalysis'];
cellname = {'GHR63T1', 'HC59T1'};
clear cellclumask;
for cl = 1: 19
    cellclumask{1, cl} = ['cluster', num2str(cl), '_mask.nii'];
end
% clusrer 1-19
clear GHRHCclusterVolnanmean;
for j= 1: length(cellname)
    clear subfiles;
    clear signal;
    subfileW = [WKD, '/', cellname{j}];
     subfiles = dir([subfileW, '/*.nii']);% image file;
    for i = 1: length(cellclumask) % read mask;
        maskvol = spm_vol([WKD, '/twosamTtestagesexeduTIV/ClusterSaveGRFvoxel0001cluster005/',...
            cellclumask{i}]);% need revise;
        maskAAL = spm_read_vols(maskvol);
        
        for subi = 1: length(subfiles)
            disp(['Processing the ', num2str(subi), '-th sub of ', cellname{j}, '!']);
            splitcell = split(subfiles(subi).name, '-');
            split1cell = split(splitcell{2}, '_');
            signal(subi, 1) = str2num(split1cell{1});% save the sub ID;
            
            V = spm_vol([subfileW, filesep, subfiles(subi).name]);
            [con1, XYZ] = spm_read_vols(V);  %  spm read image
            % con1--113*137*113double; XYZ--3*1749353double;
            con1_data = con1(maskAAL>0); %  this is a vector
            mean = nanmean(con1_data(:));
            signal(subi, i+1) = mean* 1.5^3/10^3;% mm3 
        end
    end
    GHR63HC59clusterVolmean{j, 1} = signal; 
end

% save
save  GHR63HC59cluster19VolmeanXin.mat GHR63HC59clusterVolmean;
%% 2.correcction between cluster volume and symptoms scores(SIPS/GAF); 
% using age,sex,edu,TIV as covariates;
% rowINDEX  = [1: 9, 11:63]; % Remove GHR10181;
rowINDEX  = [1:63];
% col = [27,51];% SIPS/GAF;
% col = [ 55];% LES;
col = [64, 65, 68, 72, 74, 76];
ROIco = [2: 20];% column index of cluster area vol;

clear RPGHR;
clear Cov;
Cov = GHRdataMatchImageFDTIV1{rowINDEX, [5, 4, 9, 11]};% Age,sex,edu, TIV;
for  i =1: length(ROIco) % for each cluster
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(GHR63HC59clusterVolmean{1, 1}(rowINDEX, ROIco(i)),...
            GHRdataMatchImageFDTIV1{rowINDEX, col(j)},...
            Cov, 'Type', 'Spearman', 'Rows', 'Complete');
        RPGHR(i, colu-1) = R;
        RPGHR(i, colu) = P;
    end
end

%% %%%%%%%%%%%%% GHR subgroup and HC59 %%%%%%%%%%%%%%%%%%%%%
%% extract the subgroup age,sex,edu, TIV behavior data and image data;
%% two T using age.sex.edu,TIV as COVs compare the group of GHRmax-HC59 and GHRmin-HC59;
% using GRF corrected the T statistical results and using the saved cluster
% mask of passing the GRF correction 
%% 1. to extract the cluster volume in the group of GHRmax and GHRmin, separately;
%%%%%%%%%%%%%%%  VBM 
WKD = ['~/3ndComputVBMcat12modulated/subgroup_by_HC7778vol_analysis'];
cellname = {'78RGHRmaxsubG', 'HC59T1'};
clear cellclumask;
for cl = 1: 12
    cellclumask{1, cl} = ['cluster', num2str(cl), '_mask.nii'];
end
% clusrer 1-12
clear GHRmaxHCclusterVolml;
for j= 1: length(cellname)
    clear subfiles;
    clear signal;
    subfileW = [WKD, '/', cellname{j}];
     subfiles = dir([subfileW, '/*.nii']);% image file;
    for i = 1: length(cellclumask) % read mask;
        maskvol = spm_vol([WKD, '/78GHRmax_HC59twoTagesexeduTIV/ClusterSaveGRFvoxel0001cluster005/',...
            cellclumask{i}]);% need revise;
        maskAAL = spm_read_vols(maskvol);
        
        for subi = 1: length(subfiles)
            disp(['Processing the ', num2str(subi), '-th sub of ', cellname{j}, '!']);
            splitcell = split(subfiles(subi).name, '-');
            split1cell = split(splitcell{2}, '_');
            signal(subi, 1) = str2num(split1cell{1});% save the sub ID;
            
            V = spm_vol([subfileW, filesep, subfiles(subi).name]);
            [con1, XYZ] = spm_read_vols(V);  %  spm read image
            % con1--113*137*113double; XYZ--3*1749353double;
            con1_data = con1(maskAAL>0); %  this is a vector
            mean = nanmean(con1_data(:));
            signal(subi, i+1) = mean* 1.5^3/10^3;% mm3 --> ml, %%%%%%%% mean volume;
        end
    end
    GHRmaxHCclusterVolml{j, 1} = signal; 
end
%% 2.correcction between cluster volume and symptoms scores(SIPS/GAF); 
% using age,sex,edu,TIV as covariates;
% MAX group;
rowINDEX  = GHRmaxindex78;
col = [27,51];% SIPS/GAF, sickN;
% col = [55];% LES;
ROIco = [2: 13];% column index of cluster area vol;

clear RPGHR;
clear Cov;
Cov = GHRdataMatchImageFDTIV1{rowINDEX, [5, 4, 9, 11]};% Age,sex,edu, TIV;
for  i =1: length(ROIco) % for each cluster
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(GHRmaxHCclusterVolml{1, 1}(:, ROIco(i)),...
            GHRdataMatchImageFDTIV1{rowINDEX, col(j)},...
            Cov, 'Type', 'Spearman', 'Rows', 'Complete');
        RPGHR(i, colu-1) = R;
        RPGHR(i, colu) = P;
    end
end

% MIN group;
% rowINDEX  = GHRminindex78([1:5, 7:end]);% remove 10181 for abnormal LES;
rowINDEX  = GHRminindex78;
col = [27, 51];% SIPS/GAF, sickN;
% col = [55];% LES;
ROIco = [2: 4];% column index of cluster area vol;

clear RPGHR;
clear Cov;
Cov = GHRdataMatchImageFDTIV1{rowINDEX, [5, 4, 9, 11]};% Age,sex,edu, TIV;
for  i =1: length(ROIco) % for each cluster
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(GHRminHCclusterVolml{1, 1}(:, ROIco(i)),...
            GHRdataMatchImageFDTIV1{rowINDEX, col(j)},...
            Cov, 'Type', 'Spearman', 'Rows', 'Complete');
        RPGHR(i, colu-1) = R;
        RPGHR(i, colu) = P;
    end
end

% HC59 -MINHC cluster group;
rowINDEX  = [1:11, 13: 59];
col = [26, 27];% LES;
ROIco = [2: 4];% column index of cluster area vol;

clear RPGHR;
clear Cov;
Cov = HCdataMatchImageFDTIVspm1{rowINDEX, [5, 4, 9, 11]};% Age,sex,edu, TIV;
for  i =1: length(ROIco) % for each cluster
    for j = 1:length(col)
        colu = j*2;
        [R, P] = partialcorr(GHRminHCclusterVolml{2, 1}(rowINDEX, ROIco(i)),...
            HCdataMatchImageFDTIVspm1{rowINDEX, col(j)},...
            Cov, 'Type', 'Spearman', 'Rows', 'Complete');
        RPGHR(i, colu-1) = R;
        RPGHR(i, colu) = P;
    end
end

clear fdrQ;
fdrQ = mafdr(pvalue, 'BHFDR', true);

%% %%%%%% mediation analysis %%%%
%% subgroup cluster  GHR vol-LESN-symptom(SIPS);;
% 
Rowin = GHRminindex78([1:5, 7:31, 33: end]);%%% index-subgroup;
Xa = GHRminHCclusterVolml{1, 1}(:,2);%  vol
M = zscore(Xa([1:5, 7:31, 33: end]));
X = zscore(GHRdataMatchImageFDTIV1{Rowin, 55});% LESN
Y = zscore(GHRdataMatchImageFDTIV1{Rowin, 27});% sips total
COV = GHRdataMatchImageFDTIV1{Rowin, [5, 4, 9, 11]};% age, sex, edu, TIV;

[paths, stats] = mediation(X, Y, M, 'boottop', 'bootsamples', 3000,...
    'covs', COV, 'plots', 'verbose',...
     'names', {'LESN' 'SIPStotal' 'LcaudateVol'});
 
 
 
 %% add analysis of enrichment for GHRmin and GHRmax subgroup;
 %% % read the image data Cell{sub37*1}--GHRmin-voxel mediation
 WKD = ['~/3ndComputVBMcat12modulated'];
maskvol = spm_vol([WKD, '/slice_AAL3v1nocerelluBin.nii']);
maskAAL = spm_read_vols(maskvol); % AAL Mask

    clear niifiles;
    loadpath = [WKD, filesep, '/subgroup_by_HC7778vol_analysis/78RGHRminsubG'];
    niifiles = dir([loadpath, '/*.nii']);
    for j = 1:length(niifiles)
        clear V;
        clear volume;
        V = spm_vol([loadpath, filesep, niifiles(j).name]);
        [con1, XYZ] = spm_read_vols(V);  %  spm read image
        % con1--113*137*113double; XYZ--3*1749353double;
        
        con1_data = con1(maskAAL>0); % get the voxles within the AAL mask only, this is a vector
        volumes{j, 1} = con1_data; % Cell{sub37*load1234}--382723*1double
    end

[m, n] = size(volumes); %37*1
Rowin = GHRminindex78;%%% index-subgroup;
X = GHRdataMatchImageFDTIV1{Rowin, 55};% LESN,
Y = GHRdataMatchImageFDTIV1{Rowin, 27};% sips total
COV = GHRdataMatchImageFDTIV1{Rowin, [5, 4, 9, 11]};% age, sex, edu, TIV;

for k = 1: length(volumes{1,1}) % for each location--382723
    display(['Computing the ', num2str(k), '-th location point!']);
    clear subidata;
    for row = 1: m
        subidata(row, 1) = volumes{row, 1}(k);
    end
    LOCdata{k, 1} = subidata;
    % fitlme model 1
    subidata(isnan(subidata) == 1)=[];
    Tvalue(k, 1) = 0;
    pvalueD(k, 1) = 1;

    dataROICov = [subidata(:, 1), X, Y];% volume ,LESN, sipS;
    indNan = find(sum(isnan(dataROICov), 2)>0); % find the row including nan;
   indnoNan = setdiff([1: length(Rowin)],  indNan); % 32
    dataROICov1 = zscore(dataROICov(indnoNan', :));
    
    [paths, stats] = mediation(dataROICov1(:, 2), dataROICov1(:, 3), dataROICov1(:, 1), 'boottop', 'bootsamples', 3000,...
        'covs', COV(indnoNan', :), 'verbose',...
        'names', {'LESN' 'SIPStotal' 'LcaudateVol'});
    close all;
    table1 = [stats.mean; stats.ste; stats.z; stats.ci(:, :, 1); stats.ci(:, :, 2); stats.p];
    table2 = array2table(table1, 'VariableNames', {'a', 'b', 'c1', 'c', 'ab'},...
        'RowNames', {'Coef', 'STE', 'Z', 'CI lb', 'CI ub', 'p'});
    statsSave{k, 1} = table2; % save the results;
    
    Cz(k, 1) =stats.z(4); % statistical z of variable c;
    Cp(k, 1) =stats.p(4); % statistical p of variable c(total effect);
    ABz(k, 1) =stats.z(5); % statistical z of variable ab;
    ABp(k, 1) =stats.p(5); % statistical p of variable ab(indirect path);
end
clearvars -except statsSave Cz Cp ABz ABp;
%Cz
V.fname = 'voxelMediationGHR78min_Cz.nii';
Fcon = con1;
statvalue = Cz;
Fcon(maskAAL<1)= 0;
Fcon(maskAAL>0)= statvalue;
spm_write_vol(V, Fcon);
%ABz
V.fname = 'voxelMediationGHR78min_ABz.nii';
Fcon1 = con1;
statvalue1 = ABz;
Fcon1(maskAAL<1)= 0;
Fcon1(maskAAL>0)= statvalue1;
spm_write_vol(V, Fcon1);

%%  step1:extract the mean z value in brain ROI of cerellum part of AAL3;
% in AAL3, [1:94, 121:170]; of which no label equals to 35,36,81,82;
 path = ['/public/home/Fengnana/DATA210725',...
     '/SCZ2015XYMRS/PROdcmniixBIDSnoMRS/BIDSNEWfmriprep/3ndComputVBMcat12modulated'];
mask1 = zeros(113,137,113);
Vatlas = spm_vol([path, '/slice_AAL3v1.nii']);% 1.5*1.5*1.5;
Tatlas = spm_read_vols(Vatlas); % AAL3 data;

Vmap = spm_vol([path, '/subgroup_by_HC7778vol_analysis/addAnavoxelMediationGHRmin/voxelMediationGHR78min_ABz.nii']);% 1.5*1.5*1.5;
Tmap = spm_read_vols(Vmap); % Tmap data;
label = [1:94, 121:170];
for i = 1:length(label)
    maskX = mask1+ (Tatlas == label(i));
    data = Tmap(maskX > 0);
    voxelMediationGHR78minAAL3mean(i, 1) = label(i);
    voxelMediationGHR78minAAL3mean(i, 2) = nanmean(data);
end
% save the aal3 ROI144 mean Tvalue(144*2) as csv
csvwrite('voxelMediationGHR78min_ABz_AAL3meanAAL3ROI144mean.csv', voxelMediationGHR78minAAL3mean);

%% %%%%%%%%%%%%GHRmax: similiar voxel mediation with above GHR78min subgroup%%%
%% % read the image data Cell{sub26*1}--GHRmax-voxel mediation