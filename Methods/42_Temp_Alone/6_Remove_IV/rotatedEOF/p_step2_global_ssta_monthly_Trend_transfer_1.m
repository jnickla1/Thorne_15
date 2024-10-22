load data_Global_monthly_SSTa_EOFs_1.mat              ...
        lons lats t sg cti tamo amo3 pdo       ...
        La PCa EOFa reg_EOFa stdpca eigva         

%% trend transfer algorithm
fprintf('Computing the trend transfer ...\n');

[a1,reg_a1,varex1,d1] = Trend_Transfer(PCa(:,1:4),reg_EOFa(:,:,1:4),La(1:4));

save data_Global_monthly_SSTa_trend_transfer_1.mat   ...
    a1 reg_a1 varex1 d1
%% testing code
Figures_EOFs(varex1,a1,reg_a1,lons,lats,t,[-0.5 0.5],[1910 2024],22);
