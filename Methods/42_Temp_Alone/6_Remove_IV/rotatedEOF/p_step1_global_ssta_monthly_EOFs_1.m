%% ERSSTv5: 1910-2023
fprintf('Reading ERSST v5 data ...\n');
A = load('/Volumes/chenxy/MatData/ERSSTv5/ERSST.v5.1910.2023.mat');

isyear = 1910;
ieyear = 2023;
lons = A.lons; nlons  = A.nlons;
lats = A.lats; nlats  = A.nlats;
t    = A.t;    ntimes = A.ntimes;
sg   = A.sg;
ssta = A.ssta;
ssts = A.ssts;
arw  = A.arw;
mask3 = A.mask3;

% widely used time series: cti (Deser and Wallace) and amo (amo1: Trenberth and Shea)
%
%% CTI: cti, ctis, reg_cti
fprintf('Computing CTI ...\n');
rlon2 = find(lons >= 180 & lons <= 360-90);
rlat2 = find(lats >= -6 & lats <= 6);

a1 = arw(rlon2,rlat2);
s1 = ssts(rlon2,rlat2,:);
s1 = s1.*repmat(a1,[1 1 size(s1,3)]);
cti = nansum1(reshape(s1,[size(s1,1)*size(s1,2) size(s1,3)]),1)/nansum1(a1(:));
cti = cti(:);
ctis = cti/std(cti);

reg_cti = linear_regression_2d(cti/std(cti),ssts);

%% AMO (Enfield et al. GRL, 2001): amo3, amo3s, reg_amo3
% ten-year running mean of detrended Atlantic SSTA north of Equator
fprintf('Computing AMO 3 ...\n');
sstar = detrend3d(ssta);
s1 = sstar;
s1(mask3 ~= 3) = nan;
rlats = find(lats < 0);
s1(:,rlats,:) = nan;

amo3 = area_weighted_mean(s1,lons,lats);
amo3 = runmean(amo3,13*12+1,1);
amo3s = amo3/nanstd(amo3);
rnnan = find(~isnan(amo3));
reg_amo3 = linear_regression_2d(amo3(rnnan)/std(amo3(rnnan)),ssta(:,:,rnnan));
tamo = t;

%% PDO: Manuta et al. 1997: lpdo, pdo, reg_pdo
fprintf('Computing PDO ...\n');

s1 = ssts; % Using this after the discussion with Mike.
rlat = find(lats < 20);
s1(mask3(:,:,2:end-1) ~= 2) = nan;
s1(:,rlat,:) = nan;
% weighted by sqrt of cosin of latitude
for j = 1:size(s1,2)
    s1(:,j,:) = s1(:,j,:)*sqrt(cosd(lats(j)));
end
N = 40;
[Lpdo,PCpdo,EOFpdo,stdpcpdo,eigvpdo] = myEOF3dn(s1,N);
rt = find(t >= 1997.5 & t <= 1998.5);
if mean(PCpdo(rt,1)) < 0
    PCpdo(:,1) = -PCpdo(:,1); EOFpdo(:,:,1) = -EOFpdo(:,:,1);
end
lpdo = Lpdo(1);
pdo = PCpdo(:,1);
reg_pdo = linear_regression_2d_pc(pdo,ssts);

%% Conventional EOF analysis of SSTa: La, PCa, EOFa, reg_EOFa
fprintf('Computing Conventional EOF analysis ...\n');

s1 = ssta;
% weighted by sqrt of cosin of latitude
for j = 1:size(s1,2)
    s1(:,j,:) = s1(:,j,:)*sqrt(cosd(lats(j)));
end
N = 40;
[La,PCa,EOFa,stdpca,eigva] = myEOF3dn(s1,N);
%
% check sign
% - PC1 trend, ENSO-like
% CTI region
rlon1 = find(lons >= 180 & lons <= 360-90);
rlat1 = find(lats >= -6 & lats <= 6);    
% PC1: trend
aa1 = EOFa(rlon1,rlat1,1);
aa1 = nanmean(aa1(:));
if aa1 < 0
    PCa(:,1) = -PCa(:,1);
    EOFa(:,:,1) = -EOFa(:,:,1);
end
% PC2: ENSO-LIKE
aa1 = EOFa(rlon1,rlat1,2);
aa1 = nanmean(aa1(:));
if aa1 < 0
    PCa(:,2) = -PCa(:,2);
    EOFa(:,:,2) = -EOFa(:,:,2);
end

% PC4
aa1 = EOFa(rlon1,rlat1,4);
aa1 = nanmean(aa1(:));
if aa1 > 0
    PCa(:,4) = -PCa(:,4);
    EOFa(:,:,4) = -EOFa(:,:,4);
end

% PC3, AMO-like
% North Atlantic region
rlon2 = find(lons >= 270 & lons <= 360);
rlat2 = find(lats >= 40 & lats <= 60);    
aa1 = EOFa(rlon2,rlat2,3);
aa1 = nanmean(aa1(:));
if aa1 < 0
    PCa(:,3) = -PCa(:,3);
    EOFa(:,:,3) = -EOFa(:,:,3);
end
clear reg_EOFa
for i = 1:6
    reg_EOFa(:,:,i) = linear_regression_2d(PCa(:,i),ssta);
end
%% testing code
Figures_EOFs_5(La(1:5),PCa(:,1:5),reg_EOFa(:,:,1:5),lons,lats,t,[-0.5 0.5],[1910 2024],21);

%% Output
fprintf('outputing ...\n');

save data_Global_monthly_SSTa_EOFs_1.mat            ...
        lons lats t sg cti tamo amo3 pdo       ...
        La PCa EOFa reg_EOFa stdpca eigva         
    