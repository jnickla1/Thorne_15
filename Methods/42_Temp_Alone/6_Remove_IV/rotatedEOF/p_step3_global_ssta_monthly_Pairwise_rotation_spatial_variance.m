load('data_Global_monthly_SSTa_EOFs_1.mat');
%         lons lats t sgr cti tamo amo3 pdo       ...
%         La PCa EOFa reg_EOFa stdpca eigva         

load('data_Global_monthly_SSTa_trend_transfer_1.mat');
%     a1 reg_a1 varex1 d1

%% rotation between PC3 and PC4 after the trend transfer: spatial variance
fprintf('Computing rotation based on Spatial Variance ...\n');
p2 = a1(:,3);
p4 = a1(:,4);
sp2 = reg_a1(:,:,3);
sp4 = reg_a1(:,:,4);

% Frequency rotation
[spsa,alpha1] = rotation_spatial_variance(sp2,sp4,lons,lats);

spa = spsa/std(spsa);
[minsp1,iminspa] = min(spa);
alphaspa = alpha1(iminspa);  % 35

clear aa reg_aa va
aa(:,1) = p2;
aa(:,2) = p4;
reg_aa(:,:,1) = sp2;
reg_aa(:,:,2) = sp4;
va(1) = varex1(3);
va(2) = varex1(4);
[ars1,reg_ars1,vars1] = rotation_PCs_EOFs(aa,reg_aa,va,alphaspa);

% %% testing code
% Figures_EOFs(vars1,ars1,reg_ars1,lons,lats,t,[-0.5 0.5],[1910 2017],22);

%%
fprintf('outputing 1 ...\n');
clear aa reg_aa varexa
aa(:,1) = a1(:,1);  % PC1 after trend transfer
aa(:,2) = a1(:,2);
aa(:,3) = ars1(:,2);
aa(:,4) = ars1(:,1);;

reg_aa(:,:,1) = reg_a1(:,:,1);
reg_aa(:,:,2) = reg_a1(:,:,2);
reg_aa(:,:,3) = reg_ars1(:,:,2);
reg_aa(:,:,4) = reg_ars1(:,:,1);

varexa(1) = varex1(1);
varexa(2) = varex1(2);
varexa(3) = vars1(2);
varexa(4) = vars1(1);


%% rotation between PC2 and rotated PC3 (Pacific mode): spatial variance
fprintf('Computing rotation 2 based on spatial variance ...\n');
p2 = aa(:,2);
p4 = aa(:,3);
sp2 = reg_aa(:,:,2);
sp4 = reg_aa(:,:,3);

% Frequency rotation
[spsb,alpha1] = rotation_spatial_variance(sp2,sp4,lons,lats);
spb = spsb/std(spsb);
[minsp1,iminspb] = min(spb);
alphaspb = alpha1(iminspb);  % 18
%%
clear as reg_as vas
as(:,1) = p2;
as(:,2) = p4;
reg_as(:,:,1) = sp2;
reg_as(:,:,2) = sp4;
vas(1) = varex1(3);
vas(2) = varex1(4);
[ars1,reg_ars1,vars1] = rotation_PCs_EOFs(as,reg_as,vas,alphaspb);

%% testing code
% Figures_EOFs(vars1,ars1,reg_ars1,lons,lats,t,[-0.5 0.5],[1910 2024],23);
%
%%
fprintf('outputing 2 ...\n');
clear a2 reg_a2 varex2
a2(:,1) = aa(:,1);  % PC1 after trend transfer
a2(:,2) = ars1(:,1);
a2(:,3) = ars1(:,2);
a2(:,4) = aa(:,4);

reg_a2(:,:,1) = reg_aa(:,:,1);
reg_a2(:,:,2) = reg_ars1(:,:,1);
reg_a2(:,:,3) = reg_ars1(:,:,2);
reg_a2(:,:,4) = reg_aa(:,:,4);

varex2(1) = varex1(1);
varex2(2) = vars1(1);
varex2(3) = vars1(2);
varex2(4) = varex1(4);

save data_Global_monthly_SSTa_rotation_spatial_variance_1.mat   ...
       lons lats t                                         ...
       alphaspa spa spsa iminspa                           ...       
       alphaspb spb spsb iminspb                           ...
       a2 reg_a2 varex2
   
%% rotation between PC3 and PC4 after the trend transfer: frequency cospectrum
fprintf('Computing rotation based on Frequency Cospectrum ...\n');
p2 = a1(:,3);
p4 = a1(:,4);
sp2 = reg_a1(:,:,3);
sp4 = reg_a1(:,:,4);

% Frequency rotation
[spsa,alpha1] = rotation_frequency_dependence_Lanczos(p2,p4);
spa = spsa/std(spsa);
[minsp1,iminspa] = min(spa);
alphaspa = alpha1(iminspa);  % 38

clear aa reg_aa va
aa(:,1) = p2;
aa(:,2) = p4;
reg_aa(:,:,1) = sp2;
reg_aa(:,:,2) = sp4;
va(1) = varex1(3);
va(2) = varex1(4);
[ars1,reg_ars1,vars1] = rotation_PCs_EOFs(aa,reg_aa,va,alphaspa);

% %% testing code
% Figures_EOFs(vars1,ars1,reg_ars1,lons,lats,t,[-0.5 0.5],[1910 2017],22);
%
fprintf('outputing 1 ...\n');
clear aa reg_aa varexa
aa(:,1) = a1(:,1);  % PC1 after trend transfer
aa(:,2) = a1(:,2);
aa(:,3) = ars1(:,2);
aa(:,4) = ars1(:,1);;

reg_aa(:,:,1) = reg_a1(:,:,1);
reg_aa(:,:,2) = reg_a1(:,:,2);
reg_aa(:,:,3) = reg_ars1(:,:,2);
reg_aa(:,:,4) = reg_ars1(:,:,1);

varexa(1) = varex1(1);
varexa(2) = varex1(2);
varexa(3) = vars1(2);
varexa(4) = vars1(1);


%% rotation between PC2 and rotated PC3 (Pacific mode): frequency cospectrum
fprintf('Computing rotation 2 based on Frequency Cospectrum ...\n');
p2 = aa(:,2);
p4 = aa(:,3);
sp2 = reg_aa(:,:,2);
sp4 = reg_aa(:,:,3);

% Frequency rotation
[spsb,alpha1] = rotation_frequency_dependence_Lanczos(p2,p4);
spb = spsb/std(spsb);
[minsp1,iminspb] = min(spb);
alphaspb = alpha1(iminspb);  % 49
%%
clear as reg_as vas
as(:,1) = p2;
as(:,2) = p4;
reg_as(:,:,1) = sp2;
reg_as(:,:,2) = sp4;
vas(1) = varex1(3);
vas(2) = varex1(4);
[ars1,reg_ars1,vars1] = rotation_PCs_EOFs(as,reg_as,vas,alphaspb);

% %% testing code
% Figures_EOFs(vars1,ars1,reg_ars1,lons,lats,t,[-0.5 0.5],[1910 2017],23);
%
%%
fprintf('outputing 2 ...\n');
clear a2 reg_a2 varex2
a2(:,1) = aa(:,1);  % PC1 after trend transfer
a2(:,2) = ars1(:,1);
a2(:,3) = ars1(:,2);
a2(:,4) = aa(:,4);

reg_a2(:,:,1) = reg_aa(:,:,1);
reg_a2(:,:,2) = reg_ars1(:,:,1);
reg_a2(:,:,3) = reg_ars1(:,:,2);
reg_a2(:,:,4) = reg_aa(:,:,4);

varex2(1) = varex1(1);
varex2(2) = vars1(1);
varex2(3) = vars1(2);
varex2(4) = varex1(4);

%% testing code
Figures_EOFs(varex2,a2,reg_a2,lons,lats,t,[-0.5 0.5],[1910 2024],24);
%%
save data_Global_monthly_SSTa_rotation_frequency_cospectrum_1.mat   ...
       lons lats t                                             ...
       alphaspa spa spsa iminspa                               ...       
       alphaspb spb spsb iminspb                               ...
       a2 reg_a2 varex2
