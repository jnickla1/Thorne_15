%linear CO2-temperature anomaly analysis to get HIW
%A.Jarvis 
%initiated - Jan. 2023 (following LEC380 sessions from Jan 2015 onwards)
%last edited - Sept. 2024
%no doubt this contains errors but is offered in good faith
%it would certainly benefit from some re-editing given the long and winding road that birthed it

clear all

%Law Dome ice core pCO2 data 
%https://cdiac.ess-dive.lbl.gov/ftp/trends/co2/lawdome.combined.dat
%https://www.ncei.noaa.gov/pub/data/paleo/icecore/antarctica/law/law2006.xls
%YEAR pCO2(ppm) (uncertainty quoted at +/- 1ppm undefined)
M = [2006	378.7
2005	376.7
2004	374.7
2003	372.8
2002	370.5
2001	368.3
2000	366.8
1999	365.5
1998	363.6
1997	360.3
1996	359.8
1995	358.3
1994	356.6
1993	354.9
1992	353.8
1991	352.6
1990	350.3
1989	349.5
1988	348.3
1987	346.1
1986	344.5
1985	343.2
1984	341.5
1983	340.1
1982	338.1
1981	337.6
1980	336.6
1979	334.2
1978	333.5
1977	331.7
1976	331.2
1974	328.1
1973	327.8
1972	324.1
1971	325.0
1970	324.0
1969	323.7
1967	322.9
1966	319.1
1964	318.6
1963	318.0
1961	315.7
1959	316.3
1958	314.4
1957	314.0
1956	316.3
1955	314.0
1954	311.8
1953	312.1
1950	312.6
1949	310.8
1948	310.5
1947	310.7
1946	311.5
1945	309.7
1944	311.6
1943	310.8
1942	311.3
1941	310.6
1940	311.4
1939	310.7
1938	309.6
1937	308.4
1935	308.5
1933	307.5
1929	306.7
1928	305.2
1927	305.0
1925	304.4
1924	305.2
1923	303.2
1919	303.4
1916	301.3
1914	300.4
1913	300.7
1912	298.4
1909	300.4
1906	297.7
1905	299.0
1904	295.1
1902	295.3
1901	296.5
1899	295.6
1896	298.2
1894	293.8
1893	294.6
1892	294.7
1889	291.9
1887	293.7
1886	290.6
1884	289.4
1883	291.9
1878	288.8
1874	290.5
1873	287.2
1870	287.4
1869	287.7
1867	285.2
1864	285.4
1862	286.6
1859	286.5
1855	284.9
1854	287.0
1852	288.6
1851	285.2
1849	287.7
1848	286.1
1846	284.1
1844	286.5
1841	283.0
1838	284.1
1835	283.7
1833	284.5
1827	285.1
1826	281.3
1814	284.3
1800	281.1
1799	283.7
1796	281.6
1794	281.5
1781	276.8
1780	279.5
1774	277.8
1763	276.5
1752	276.8
1749	276.9
1743	276.7
1734	278.2
1723	277.2
1694	276.5
1690	276.3
1682	275.9
1649	277.2
1640	276.6
1629	274.5
1610	271.8
1603	274.3
1591	278.7
1588	281.0
1573	281.9
1560	281.7
1550	282.8
1530	283.2
1502	282.4
1469	279.6
1449	281.7
1431	282.5
1429	279.5
1411	279.6
1391	280.0
1390	280.4
1350	280.1
1330	283.4
1306	281.5
1276	281.1
1258	282.1
1246	281.7
1207	283.6
1193	283.9
1160	283.9
1137	283.8
1105	282.8
1088	282.4
1058	282.8
1037	280.3
1025	280.8
1005	279.9
968	278.5
944	279.1
897	278.9
857	279.3
799	278.5
764	278.5
730	278.5
698	279.7
668	279.4
632	278.3
596	276.9
572	277.6
537	276.0
500	276.4
461	276.7
428	276.9
365	277.0
329	278.9
302	279.8
274	280.1
228	281.5
202	280.7
168	280.1
136	278.1
104	277.5
56	277.4
30	277.9
13	276.7];
yearLawDome = flipud(round(M(:,1)));
LawDome = flipud(M(:,2));

%Mauna Loa atmospheric pCO2 data
%https://gml.noaa.gov/webdata/ccgg/trends/co2/co2_annmean_mlo.txt
%YEAR [pCO2](ppm) 'UNCERTAINTY'(UNDEFINED?)
M = [1959   315.98     0.12
  1960   316.91     0.12
  1961   317.64     0.12
  1962   318.45     0.12
  1963   318.99     0.12
  1964   319.62     0.12
  1965   320.04     0.12
  1966   321.37     0.12
  1967   322.18     0.12
  1968   323.05     0.12
  1969   324.62     0.12
  1970   325.68     0.12
  1971   326.32     0.12
  1972   327.46     0.12
  1973   329.68     0.12
  1974   330.19     0.12
  1975   331.13     0.12
  1976   332.03     0.12
  1977   333.84     0.12
  1978   335.41     0.12
  1979   336.84     0.12
  1980   338.76     0.12
  1981   340.12     0.12
  1982   341.48     0.12
  1983   343.15     0.12
  1984   344.87     0.12
  1985   346.35     0.12
  1986   347.61     0.12
  1987   349.31     0.12
  1988   351.69     0.12
  1989   353.20     0.12
  1990   354.45     0.12
  1991   355.70     0.12
  1992   356.54     0.12
  1993   357.21     0.12
  1994   358.96     0.12
  1995   360.97     0.12
  1996   362.74     0.12
  1997   363.88     0.12
  1998   366.84     0.12
  1999   368.54     0.12
  2000   369.71     0.12
  2001   371.32     0.12
  2002   373.45     0.12
  2003   375.98     0.12
  2004   377.70     0.12
  2005   379.98     0.12
  2006   382.09     0.12
  2007   384.02     0.12
  2008   385.83     0.12
  2009   387.64     0.12
  2010   390.10     0.12
  2011   391.85     0.12
  2012   394.06     0.12
  2013   396.74     0.12
  2014   398.81     0.12
  2015   401.01     0.12
  2016   404.41     0.12
  2017   406.76     0.12
  2018   408.72     0.12
  2019   411.65     0.12
  2020   414.21     0.12
  2021   416.41     0.12
  2022   418.53     0.12
  2023   421.08     0.12];
yearMaunaLoa = M(:,1);
MaunaLoa = M(:,2);

%HadCRUT5A global temperature anomaly 
%https://www.metoffice.gov.uk/hadobs/hadcrut5/data/current/analysis/diagnostics/HadCRUT.5.0.1.0.analysis.summary_series.global.annual.csv
%YEAR HadCRUT5A(degC) 5thPERCENTILE(degC) 95thPERCENTILE(degC)
M = [1850	-0.4177114	-0.58925647	-0.2461663
1851	-0.2333498	-0.41186792	-0.054831687
1852	-0.22939907	-0.40938243	-0.04941572
1853	-0.27035445	-0.43000934	-0.110699534
1854	-0.29152083	-0.4327115	-0.15033019
1855	-0.29691675	-0.43931878	-0.15451474
1856	-0.32035372	-0.46809322	-0.1726142
1857	-0.46723005	-0.61632216	-0.31813794
1858	-0.3887657	-0.53688604	-0.24064532
1859	-0.28126517	-0.42392084	-0.1386095
1860	-0.39016518	-0.5389766	-0.24135375
1861	-0.42911294	-0.597001	-0.26122484
1862	-0.5363694	-0.70365995	-0.3690788
1863	-0.34424406	-0.53403676	-0.15445141
1864	-0.46546507	-0.64809304	-0.28283712
1865	-0.33248132	-0.5244569	-0.14050578
1866	-0.3412875	-0.5218009	-0.16077414
1867	-0.35699412	-0.5530232	-0.16096504
1868	-0.35182714	-0.5294787	-0.1741756
1869	-0.31659195	-0.47646335	-0.15672056
1870	-0.32792753	-0.46869946	-0.18715557
1871	-0.36856276	-0.51413316	-0.22299235
1872	-0.32811058	-0.46317002	-0.19305114
1873	-0.3412969	-0.47250164	-0.21009217
1874	-0.3732512	-0.5071426	-0.2393598
1875	-0.37562594	-0.514041	-0.23721085
1876	-0.4241099	-0.56287116	-0.28534868
1877	-0.10110884	-0.22982001	0.02760234
1878	-0.011315192	-0.13121258	0.1085822
1879	-0.30363432	-0.43406436	-0.1732043
1880	-0.31583205	-0.44015092	-0.19151321
1881	-0.23224552	-0.35793498	-0.10655606
1882	-0.29553008	-0.4201501	-0.17091006
1883	-0.3464744	-0.4608177	-0.23213112
1884	-0.49232006	-0.6026686	-0.38197154
1885	-0.47112358	-0.5830682	-0.35917896
1886	-0.42090362	-0.5225382	-0.31926903
1887	-0.49878576	-0.61655986	-0.38101166
1888	-0.37937889	-0.49332377	-0.265434
1889	-0.24989556	-0.37222093	-0.12757017
1890	-0.50685817	-0.6324095	-0.3813068
1891	-0.40131494	-0.5373699	-0.26525998
1892	-0.5075585	-0.64432853	-0.3707885
1893	-0.49461922	-0.6315314	-0.35770702
1894	-0.48376393	-0.6255681	-0.34195974
1895	-0.4487516	-0.58202064	-0.3154826
1896	-0.28400728	-0.41740146	-0.15061307
1897	-0.25980017	-0.39852425	-0.12107606
1898	-0.48579213	-0.6176492	-0.35393503
1899	-0.35543364	-0.48639694	-0.22447036
1900	-0.2344939	-0.3664699	-0.1025179
1901	-0.29341024	-0.42826617	-0.15855433
1902	-0.43895653	-0.57330817	-0.3046049
1903	-0.5332871	-0.65805364	-0.40852055
1904	-0.59751105	-0.72495615	-0.47006598
1905	-0.40779322	-0.53042954	-0.2851569
1906	-0.31910878	-0.4449256	-0.19329198
1907	-0.5040763	-0.61996025	-0.38819233
1908	-0.5138197	-0.6301149	-0.39752448
1909	-0.53568715	-0.64385575	-0.42751858
1910	-0.5309095	-0.6454714	-0.4163476
1911	-0.539079	-0.6510565	-0.42710155
1912	-0.47553864	-0.5783255	-0.37275177
1913	-0.4670111	-0.5772189	-0.3568033
1914	-0.26243657	-0.37001243	-0.15486072
1915	-0.19167219	-0.30920455	-0.07413985
1916	-0.42002314	-0.5453339	-0.2947124
1917	-0.5428197	-0.67724293	-0.40839645
1918	-0.4243641	-0.56627965	-0.28244856
1919	-0.32528907	-0.46366385	-0.18691428
1920	-0.29835507	-0.4310807	-0.16562946
1921	-0.24044435	-0.36450884	-0.11637985
1922	-0.3390137	-0.449683	-0.22834438
1923	-0.31768188	-0.4260641	-0.20929964
1924	-0.3118017	-0.4204113	-0.2031921
1925	-0.28214198	-0.39495143	-0.16933253
1926	-0.122555	-0.23084709	-0.014262918
1927	-0.2291136	-0.33219776	-0.12602946
1928	-0.20646581	-0.31743687	-0.09549474
1929	-0.39244303	-0.50286	-0.28202602
1930	-0.17680542	-0.29041144	-0.06319938
1931	-0.103397675	-0.2126916	0.005896252
1932	-0.14546171	-0.25195518	-0.03896824
1933	-0.32234442	-0.4271004	-0.21758842
1934	-0.17433685	-0.27400395	-0.07466974
1935	-0.20605923	-0.30349734	-0.10862113
1936	-0.16952094	-0.26351926	-0.07552263
1937	-0.019198947	-0.11975877	0.08136088
1938	-0.012200737	-0.110303745	0.08590227
1939	-0.040797204	-0.14670469	0.06511029
1940	0.07593582	-0.041949686	0.19382133
1941	0.0381293	-0.16225392	0.23851252
1942	0.001406068	-0.19521242	0.19802456
1943	0.006421582	-0.19958311	0.21242628
1944	0.14410514	-0.05449483	0.3427051
1945	0.043088354	-0.1572829	0.24345961
1946	-0.11881461	-0.26595104	0.028321814
1947	-0.09120561	-0.23179047	0.049379256
1948	-0.12466127	-0.25913337	0.009810847
1949	-0.14380223	-0.25407746	-0.033526964
1950	-0.2266218	-0.33265698	-0.120586626
1951	-0.06115396	-0.15035023	0.028042309
1952	0.015354548	-0.08293598	0.11364508
1953	0.07763075	-0.020529611	0.17579111
1954	-0.11675023	-0.20850272	-0.02499773
1955	-0.19730994	-0.28442997	-0.11018991
1956	-0.26316562	-0.33912566	-0.18720558
1957	-0.035334915	-0.10056861	0.029898778
1958	-0.017632563	-0.08307456	0.04780944
1959	-0.048004813	-0.11036374	0.014354111
1960	-0.11545958	-0.17413756	-0.056781605
1961	-0.019999769	-0.07077983	0.030780297
1962	-0.06404272	-0.117309235	-0.010776198
1963	-0.036810614	-0.09057977	0.016958548
1964	-0.30586153	-0.34949452	-0.26222858
1965	-0.20442048	-0.2536107	-0.15523025
1966	-0.1488976	-0.19840918	-0.09938603
1967	-0.117539294	-0.1606568	-0.07442179
1968	-0.16864756	-0.21327549	-0.12401963
1969	-0.03138624	-0.071891375	0.009118891
1970	-0.08506408	-0.12605265	-0.044075526
1971	-0.20588905	-0.24447224	-0.16730586
1972	-0.09379131	-0.1316835	-0.055899136
1973	0.04995016	0.013479018	0.0864213
1974	-0.17252657	-0.21022564	-0.13482748
1975	-0.110754214	-0.15130511	-0.07020333
1976	-0.2158369	-0.25586063	-0.17581315
1977	0.1030885	0.06005668	0.14612031
1978	0.005255972	-0.034596615	0.045108557
1979	0.09085814	0.062358625	0.11935765
1980	0.19607204	0.16280398	0.2293401
1981	0.25001204	0.21939126	0.28063282
1982	0.034268282	-0.005108635	0.0736452
1983	0.22380984	0.18804431	0.25957537
1984	0.04799352	0.011540157	0.08444688
1985	0.049729742	0.015663432	0.08379605
1986	0.09568698	0.064408004	0.12696595
1987	0.2430264	0.21218552	0.27386728
1988	0.2821517	0.24703526	0.31726816
1989	0.1792503	0.14449841	0.21400218
1990	0.36058238	0.32456872	0.39659604
1991	0.33889654	0.30403617	0.3737569
1992	0.12489683	0.09088211	0.15891157
1993	0.16570719	0.12821673	0.20319764
1994	0.23354979	0.19841295	0.26868662
1995	0.37686613	0.34365577	0.41007653
1996	0.27668938	0.24318002	0.31019875
1997	0.4223085	0.39009085	0.4545262
1998	0.5773417	0.54306877	0.61161464
1999	0.32448497	0.29283476	0.3561352
2000	0.33108476	0.29822785	0.36394167
2001	0.48928034	0.4580683	0.5204924
2002	0.5434665	0.51278186	0.57415116
2003	0.54417014	0.5112426	0.5770977
2004	0.46737072	0.43433833	0.5004031
2005	0.6068625	0.5757053	0.63801974
2006	0.5725527	0.541973	0.60313237
2007	0.5917013	0.56135315	0.62204945
2008	0.46564984	0.43265733	0.49864236
2009	0.5967816	0.5652556	0.62830764
2010	0.6803714	0.6490759	0.7116668
2011	0.53769773	0.5060012	0.5693943
2012	0.57760704	0.54485524	0.61035883
2013	0.6235753	0.5884838	0.6586669
2014	0.67287165	0.63890487	0.7068384
2015	0.8251144	0.7912871	0.85894173
2016	0.9329271	0.9017635	0.96409065
2017	0.84517425	0.81477475	0.87557375
2018	0.76265407	0.73105204	0.79425603
2019	0.8910726	0.85678726	0.9253579
2020	0.9229205	0.8882981	0.95754296
2021	0.76190555	0.7254577	0.7983534
2022	0.8013053	0.76525766	0.837353
2023	1.0858556	1.0421426	1.1295687];
yearHadCRUT5A = M(:,1);
HadCRUT5A = M(:,2);
HadCRUT5AL = M(:,3);
HadCRUT5AU = M(:,4);

%temp measurement uncertainties
HadCRUT5Arang = HadCRUT5AU-HadCRUT5AL;
HadCRUT5Asig = HadCRUT5Arang/4;  %95th percentile hence assuming +/- 2 sig. but need to explore those HadCRUT5 ensembles

%1850-1900 baseline
tempanomalybaseline = mean(HadCRUT5A(1:50));
tempanomalybaselinesig = std(HadCRUT5A(1:50)) + mean(HadCRUT5Asig(1:50));

%-------------------------------
%ANALYSIS
%-------------------------------
%construct interpolated single pCO2 record
yearpCO2 = (min(yearLawDome):max(yearMaunaLoa))';
pCO2 = zeros(length(yearpCO2),1)*NaN;

%fill in Law Dome data
for i = 1:length(yearpCO2)
    
    j = find(yearLawDome == yearpCO2(i));
    if isempty(j)==1
        pCO2(i,1) = NaN;
    elseif isempty(j)==0
        pCO2(i) = LawDome(j);
    end
    
end

%fill in Mauna Loa data
i = find(yearpCO2>=min(yearMaunaLoa));
pCO2(i) = MaunaLoa;

%how many points are interpolated?
i = find(yearpCO2>=1850 & yearpCO2<=1958);
j = find(isnan(pCO2(i))==1);
R = length(j)/length(i);

%interpolate missing pCO2 data and smooth ice core data
i = find(isnan(pCO2)==0);
[pCO2i] = interp1(yearpCO2(i),pCO2(i),yearpCO2,'pchip');%plot(yearpCO2,pCO2,'ko',yearpCO2,pCO2i,'r');return
[nvr,opts,parse] = irwsmopt([pCO2i],1,'ml');
[pCO2si] = irwsm(pCO2i,1,nvr/1e3);  %approximates Enting's 20 year spline smoothing

%retain pre 1958 data
i = find(isnan(pCO2)==0);
pCO2si(i) = pCO2(i);

%retain post 1958 data
i = find(yearpCO2>1958);
pCO2si(i) = pCO2(i);

%CO2 measurement uncertainties
pCO2sig = zeros(size(pCO2));
i = find(yearpCO2<1959);
pCO2sig(i) = 2/2;          %assume Etheridge quoted uncertainty of +/- 2ppmv - assume +/- 2sig?
i = find(yearpCO2>=1959);
pCO2sig(i) = 0.12/2;        %assume Keeling quoted uncertainty of +/- 0.12ppmv - assume +/- 2sig?

%observed CO2 baselines
i = find(yearLawDome<=1700);
pCO2baselinemu = mean(LawDome(i));
pCO2baselinesig = std(LawDome(i)) + 1;   %add Etheridge quoted measurement uncertainty of +/- 2ppmv - assume +/- 2sig   
disp('pre-1700 CO2 baseline')
disp([pCO2baselinemu pCO2baselinesig*2])

i = find(yearpCO2>=1850 & yearpCO2<=1900);
pCO2baselinemu1850_1900 = mean(pCO2si(i));
nvr = irwsmopt(pCO2si(i),1);
pCO2baselinetrend = irwsm(pCO2si(i),1,nvr);
pCO2baselinesig1850_1900 = std(pCO2si(i)) + 1 - std(pCO2baselinetrend);  %add Etheridge quoted measurement uncertainty of +/- 2ppmv - assume +/- 2sig but subtract trend variance   
disp('1850-1900 CO2 baseline')
disp([pCO2baselinemu1850_1900 pCO2baselinesig1850_1900*2])

%difference between the pre-1700 and 1850-1900 CO2 baselines
N = 1e4;
R1 = randn(N,1)*pCO2baselinesig+pCO2baselinemu;
R2 = randn(N,1)*pCO2baselinesig1850_1900+pCO2baselinemu1850_1900;
CO2baselinediffmu = pCO2baselinemu1850_1900-pCO2baselinemu;
CO2baselinediffsig = std(R2-R1);

%print baselines
disp('difference between1850-1900 CO2 baseline')
disp([CO2baselinediffmu CO2baselinediffsig*2])
disp('1850-1900 HadCRUT5 baseline')
disp([tempanomalybaseline tempanomalybaselinesig*2])

%----------------------------
% REGRESSION 1 - 1959-present
%1. define x, y and weights
i = find(yearpCO2>=1959);
x = pCO2si(i)-pCO2baselinemu;
i = find(yearHadCRUT5A>=1959);
y = HadCRUT5A(i);
wy = 1./(HadCRUT5Asig(i)).^2; 
%2. iterative WLS
A = [1 0];
for i = 1:4
    xf = filter(A,1,x);
    yf = filter(A,1,y);
    [Xf,STDX,MSE,Sf] = lscov([xf ones(size(x))],yf,wy);
    m = Xf(1);
    c = Xf(2)/sum(A);
    Xfsig = sqrt(diag(Sf));
    msig = Xfsig(1);
    csig = Xfsig(2)/sum(A);
    yhat = m*x + c;
    e = y-yhat;
    %th = ar(e,1);A = polydata(th);
    th = mar(e,1);
    [A,B,C,P] = getpar(th);
end
disp('post 1958 regression result')
disp([[m msig*2]*100;c csig*2])

% REGRESSION 2 - 1850-1958 data
%1. define x, y and weights
i = find(yearpCO2>=1850 & yearpCO2<=1958);
x = pCO2si(i)-pCO2baselinemu;   
i = find(yearHadCRUT5A>=1850 & yearHadCRUT5A<=1958);
y = HadCRUT5A(i);
wy = 1./(HadCRUT5Asig(i)).^2; 
%2. iterative WLS
A = [1 0];
for i = 1:4
    xf = filter(A,1,x);
    yf = filter(A,1,y);
    [Xf,STDX,MSE,Sf] = lscov([xf ones(size(x))],yf,wy);
    m = Xf(1);
    c = Xf(2)/sum(A);
    Xfsig = sqrt(diag(Sf));
    msig = Xfsig(1);
    csig = Xfsig(2)/sum(A);
    yhat = m*x + c;
    e = y-yhat;
    %th = ar(e,1);A = polydata(th);
    th = mar(e,1);
    [A,B,C,P] = getpar(th);
end
disp('pre 1958 regression result')
disp([[m msig*2]*100;c csig*2])

% REGRESSION 3 - 1850 - 2023
%1. define x, y and weights
i = find(yearpCO2>=1850);
x = pCO2si(i)-pCO2baselinemu;   %we could feed the baseline uncertainty into the regression here
i = find(yearHadCRUT5A>=1850);
y = HadCRUT5A(i);
wy = 1./(HadCRUT5Asig(i)).^2;
%2. iterative WLS
A = [1 0];
for i = 1:4
    xf = filter(A,1,x);
    yf = filter(A,1,y);
    [Xf,STDX,MSE,Sf] = lscov([xf ones(size(x))],yf,wy);
    m = Xf(1);
    c = Xf(2)/sum(A);
    Xfsig = sqrt(diag(Sf));
    msig = Xfsig(1);
    csig = Xfsig(2)/sum(A);
    yhat = m*x + c;
    e = y-yhat;
    %th = ar(e,1);A = polydata(th);
    th = mar(e,1);
    [A,B,C,P] = getpar(th);
end
AA = A;
PP = P;
disp('1850-2023 regression result')
disp([[m msig*2]*100;c csig*2])
AAsig = PP^0.5;
disp('AR(1) result')
disp([AA AAsig*2])

%finalised 1850-2023 estimates of senstivity and offset
MM = m;MMsig = msig;
CC = c;CCsig = csig;

%test residuals for normality
ef = filter(A,1,e); %decorrelate
efsig = std(ef);
[H,P] = adtest(ef,'alpha',0.001);
disp('adtest result')
disp([H P])

%test residuals for unit root
[h,pValue,stat,cValue,reg] = adftest(e);
disp('e is not I(1)-which you know given ar(1)<<1')

%save data for offline analysis
save CO2tempdata x y wy MM CC MMsig CCsig AA AAsig PP efsig pCO2baselinesig

%------------------------------
%Recursive WLS to find when the regression estimates stabilise?
%Produces data for Fig.2
M = 123; %start from 1900 as pre-1900 regression unstable
RECURSIVERESULT = zeros(M,7);

for j = 1:M

%1. define x, y and weights
i = find(yearpCO2==2023-M+j-1);endyear = yearpCO2(i);
i = find(yearpCO2>=1850 & yearpCO2<=endyear);
x = pCO2si(i)-pCO2baselinemu;   
i = find(yearHadCRUT5A>=1850 & yearHadCRUT5A<=endyear);
y = HadCRUT5A(i);
wy = 1./(HadCRUT5Asig(i)).^2;
%2. iterative WLS
A = [1 0];
for i = 1:4
    xf = filter(A,1,x);
    yf = filter(A,1,y);
    [Xf,STDX,MSE,Sf] = lscov([xf ones(size(x))],yf,wy);
    m = Xf(1);
    c = Xf(2)/sum(A);
    Xfsig = sqrt(diag(Sf));
    msig = Xfsig(1);
    csig = Xfsig(2)/sum(A);
    yhat = m*x + c;
    e = y-yhat;
    %th = ar(e,1);A = polydata(th);
    th = mar(e,1);
    [A,B,C,P] = getpar(th);
end

RECURSIVERESULT(j,:) = [endyear m msig c csig A(2) P^0.5];

end

%-----------------------------
%numerical baseline comparison
N = 1e4;
%R1 = HadCRUT5A(1:50);%randn(N,1)*tempanomalybaselinesig+tempanomalybaseline;
R1 = randn(N,1)*tempanomalybaselinesig+tempanomalybaseline;
R2 = randn(N,1)*CCsig + CC;
[H,P] = ttest(R1,R2);
disp('GMSTC baseline comparison')
disp([tempanomalybaseline tempanomalybaselinesig*2;...
    CC CCsig*2])
disp([CC-tempanomalybaseline std(R2-R1)*2])
disp(P)

%-----------------------------
%TCRE comparison
%should really include AF and RF uncertainties, but argument is somewhat circular
TCRE = 1e3*0.8*0.44*MM/2.123;
TCREsig = 1e3*0.8*0.44*MMsig/2.123;
disp('TCRE')
disp([TCRE TCREsig*2])

%-----------------------------
%simulate full stochastic model explicitly as don't know the analytical
%solution for CIs using weights AND stochastic baseline
%adding temp. weight effects on post MCS
%start from 1700
i = find(yearpCO2>=1700);
M = length(i);
N = 1e4;
yearMCS = yearpCO2(i);
pCO2MCS = pCO2si(i);
pCO2sigMCS = pCO2sig(i);
HadCRUT5AMCS = [ones(length(pCO2MCS)-length(HadCRUT5A),1)*NaN;HadCRUT5A];
HadCRUT5AsigMCS = [ones(length(pCO2MCS)-length(HadCRUT5A),1)*mean(HadCRUT5Asig(1:20));HadCRUT5Asig];
%deterministic x,y
HICO2D = pCO2MCS-pCO2baselinemu;
GMSTCD = MM*HICO2D;
GMSTCDO = [ones(length(pCO2MCS)-length(HadCRUT5A),1)*NaN;HadCRUT5A-CC];
HICO2MCS = zeros(M,N);
eMCS = zeros(M,N);
ECMCS = zeros(M,1);
ETMCS = zeros(M,1);
GMSTAMCS = zeros(M,N);
GMSTCMCS = zeros(M,N);
GMSTC2MCS = zeros(M,N);
SSMCS = zeros(M,N);
SOMCS = zeros(M,N);
HIWMCS = zeros(M,N);

%sample parameter covariance
T1 = chol(Sf);
N1 = randn(2,N);
pse = (T1'*N1)';
aMCS = AA(2)+randn(N,1)*AAsig;
mMCS = MM + pse(:,1);
cMCS = CC + pse(:,2)./(1+aMCS); %could include AR(1) uncertainty

%MCS iterate
for i = 1:N
    
    %check if CO2 baseline uncertainty active
    HICO2MCS(:,i) = pCO2MCS - pCO2baselinemu + randn*pCO2baselinesig*1 + randn(M,1).*pCO2sigMCS;
    eMCS(:,i) = filter(1,[1 aMCS(i)],randn(M,1)*efsig).*HadCRUT5AsigMCS./nanmean(HadCRUT5AsigMCS);   %AR(1) stochasticity weighted by HadCRUT5 uncertainty
    GMSTAMCS(:,i) = mMCS(i)* HICO2MCS(:,i) + cMCS(i) + eMCS(:,i);                   %full linear stochatic temperature anomolies
    GMSTCMCS(:,i) = GMSTAMCS(:,i) - CC;                                             %full linear stochastic GMST
    HIWMCS(:,i) =  mMCS(i)* HICO2MCS(:,i) + cMCS(i) - CC;                           %HIW trend
    SOMCS(:,i) = (GMSTCDO+GMSTCDO-GMSTCMCS(:,i))./HICO2MCS(:,i);

end

%mu's (parametric)
HICO2mu = pCO2MCS - pCO2baselinemu;
HIWmu = MM*HICO2mu;

%std's
HICO2sig = std(HICO2MCS')';
HIWsig = std(HIWMCS')';
GMSTCsig = std(GMSTCMCS')';

%percentiles
%for HIW and HICO2 use parametric mean +/-2 sig for consistency with regression results
%otherwise use non-parametric percentiles through sorting
%differences are otherwise small
HICO205 = HICO2mu-2*HICO2sig;
HICO295 = HICO2mu+2*HICO2sig;

GMSTA05 = prctile(GMSTAMCS',05)';
GMSTA95 = prctile(GMSTAMCS',95)';

GMSTC05 = prctile(GMSTCMCS',05)';
GMSTC95 = prctile(GMSTCMCS',95)';

HIW05 = HIWmu - 2*HIWsig;
HIW95 = HIWmu + 2*HIWsig;

SO05 = prctile(SOMCS',05)';
SO95 = prctile(SOMCS',95)';

%HIW 1850 - 1900
i = find(yearMCS>=1850 & yearMCS<=1900);
HIW1850_1900trendsig = std(mean(HIWMCS(i,:)')');    %remove trend variance effect for consistency
HIW1850_1900mu = mean(reshape(HIWMCS(i,:),N*length(i),1));
HIW1850_1900sig = std(reshape(HIWMCS(i,:),N*length(i),1))-HIW1850_1900trendsig;
disp('HIW 1850 - 1900')
disp([HIW1850_1900mu HIW1850_1900sig*2])

%HICO2 2023
HICO22023mu = pCO2(end)-pCO2baselinemu;
HICO22023sig = HICO2sig(end);
disp('2023 HICO2')
disp([HICO22023mu 2*HICO22023sig])

%HIW 2023
HIW2023mu = MM*(pCO2si(end) - pCO2baselinemu);  %i.e. the parametric estimate for consistency
HIW2023sig = HIWsig(end);                %+/- 2 sigma is not quite the same as the 5-95th percentile range from the MCS so Fig.1a and d do not perfectly align
HIW2023Bmu = HIW2023mu - HIW1850_1900mu;        %2023 HIW estimate baselined against 1850-1900 
HIW2023Bsig = 0.038;                           %this is the uncertainty for the 2023 HIW if we don't include the CO2 baseline uncertainty in the MCS
disp('2023 HIW')
disp([HIW2023mu 2*HIW2023sig])
disp('2023 HIW using 1850-1900 as baseline')
disp([HIW2023Bmu 2*HIW2023Bsig])

%probability that 1.5 HIW is exceeded
[HIW2023s,i]=sort(HIWMCS(end,:));
p = (1:length(i))'./(1+length(i));
%figure(3);clf;plot(GMSTC2023s,p);grid
i = find(HIW2023s>=1.5);
P15 = 1-p(i(1));
disp('P>1.5 degC')
disp(P15)

%--------------------------------------------------------------------------
%PLOTS
%--------------------------------------------------------------------------
%COL = [186 124 124]/255;
COL = [175 100 100]/255;

%-----------------------------
%Fig.1b - CO2 time series plot
mm = 50/1850;cc = 1800;    %to remap the 0AD-1700AD data
xmin = 1800;xmax = 2030;xinc = 250;
ymin = 250;ymax = 450;yinc = 50;
figure(1);clf
hold on
plot(yearMaunaLoa,MaunaLoa,'k+','markersize',4)
i = find(yearLawDome<=1958 & yearLawDome>=1850);
plot(yearLawDome(i),LawDome(i),'k.','markersize',7)
i = find(yearLawDome<=1850);
plot(yearLawDome(i)*mm+cc,LawDome(i),'k.')
plot([1850 1850],[ymin ymax],'k--')
plot([1700 1700]*mm+cc,[ymin ymax],'k--')
i = find(yearpCO2>=1850);
h = patch([yearpCO2(i);flipud(yearpCO2(i))],[pCO2si(i)+pCO2sig(i)*2;flipud(pCO2si(i)-pCO2sig(i)*2)],[1 1 1]);set(h,'facecolor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'facealpha',0.35,'edgealpha',0.35);
i = find(yearpCO2<1850);
h = patch([yearpCO2(i);flipud(yearpCO2(i))]*mm+cc,[pCO2si(i)+pCO2sig(i)*2;flipud(pCO2si(i)-pCO2sig(i)*2)],[1 1 1]);set(h,'facecolor',[0.5 0.5 0.5],'edgecolor',[0.5 0.5 0.5],'facealpha',0.35,'edgealpha',0.35);
hold off
flab = gca;
flab.XMinorTick = 'on';
flab.YMinorTick = 'on';
xtickformat('%.0f');
ytickformat('%.0f');
flab.XTick = [1850:50:2000];
flab.XTickLabel = [1850:50:2000];
%flab.YTick = [1 10 100 1000 10000];
%flab.YTickLabel = [10 10 10 10 10];
flab.FontSize = 30;
flab.FontName = 'helvetica';
flab.FontWeight = 'normal';
%xlabel('year','fontsize',40,'fontname','helvetica')
%ylabel('atmospheric [CO ] (ppmv)','fontsize',40,'fontname','helvetica')
box('on')
axis([xmin xmax ymin ymax])
%axis('square')

%--------------------------------------------------------------------------
%Fig.1a - Anomaly/GMSTc/HIW time series plot
xmin = 1800;xmax = 2030;xinc = 250;
ymin = -1;ymax = 2;yinc = 0.2;
i = find(yearMCS>=1850);
figure(2);clf
hold on
%GTA
plot(yearHadCRUT5A,HadCRUT5A,'k')
h = patch([yearHadCRUT5A;flipud(yearHadCRUT5A)],[HadCRUT5AL;flipud(HadCRUT5AU)],[1 1 1]);set(h,'facecolor',[1 1 1]*.5,'edgecolor',[1 1 1]*.5,'facealpha',0.5,'edgealpha',0.5);
%GMSTc
plot(yearMCS(i),GMSTCDO(i),'color',COL*.8)
h = patch([yearMCS(i);flipud(yearMCS(i))],[(GMSTCDO(i)-2*GMSTCsig(i));flipud(GMSTCDO(i)+2*GMSTCsig(i))],[1 1 1]);set(h,'facecolor',COL,'edgecolor',COL,'facealpha',0.35,'edgealpha',0.35);
%h = patch([yearMCS(i);flipud(yearMCS(i))],[GMSTC205(i);flipud(GMSTC295(i))],[1 1 1]);set(h,'facecolor',COL,'edgecolor',COL,'facealpha',0.25,'edgealpha',0.25);
%plot(yearMCS(i),GMSTC205(i),'k:')
%plot(yearMCS(i),GMSTC295(i),'k:')
%HIW
h = patch([yearMCS(i);flipud(yearMCS(i))],[HIW05(i);flipud(HIW95(i))],[1 1 1]);set(h,'facecolor',COL,'edgecolor',COL,'facealpha',0.55,'edgealpha',0.55);
plot([1850 1850],[ymin ymax],'k--')
plot([1900 1900],[ymin ymax],'k--')
plot([xmin xmax],[0 0],'k--')
hold off
flab = gca;
flab.XMinorTick = 'on';
flab.YMinorTick = 'on';
xtickformat('%.0f');
ytickformat('%.1f');
flab.XTick = [1850:50:2000];
flab.XTickLabel = '';%[1850:50:2000];
%flab.YTick = [1 10 100 1000 10000];
%flab.YTickLabel = [10 10 10 10 10];
flab.FontSize = 30;
flab.FontName = 'helvetica';
flab.FontWeight = 'normal';
%xlabel('year','fontsize',40,'fontname','helvetica')
%ylabel('atmospheric [CO ] (ppmv)','fontsize',40,'fontname','helvetica')
box('on')
axis([xmin xmax ymin ymax])
%axis('square')

%--------------------------------------------------------------------------
%Fig.1c - CO2-temperature plot
xmin = -50;xmax = 200;xinc = 50;
ymin = -1;ymax = 2;yinc = 0.5;
figure(3);clf
hold on
i = find(yearMCS<=1958 & yearMCS>=1850);
j = find(yearHadCRUT5A<=1958);
plot(HICO2D(i),HadCRUT5A(j),'k.','markersize',10)
i = find(yearMCS>1958);
j = find(yearHadCRUT5A>1958);
plot(HICO2D(i),HadCRUT5A(j),'k+','markersize',4)
i = find(yearMCS>=1849);
h = patch([HICO2mu(i);flipud(HICO2mu(i))],[GMSTA05(i);flipud(GMSTA95(i))],[1 1 1]);set(h,'facecolor',[1 1 1]*.5,'edgecolor',[1 1 1]*.5,'facealpha',0.35,'edgealpha',0.35);
i = find(HICO2D>=0);
plot([HICO2D(i)],[MM*HICO2D(i)+CC],'color',[0 0 0],'linewidth',1)
i = find(yearMCS>=1849);
h = patch([HICO2mu(i);flipud(HICO2mu(i))],[HIW05(i);flipud(HIW95(i))],[1 1 1]);set(h,'facecolor',COL,'edgecolor',COL,'facealpha',0.75,'edgealpha',0.75);
%plot([0;HICO2(i)],[0;MM*HICO2(i)],'color',[0.5 0 0],'linewidth',1)
plot([xmin 0],[CC CC],'k--')
plot([xmin 0],[0 0],'k--')
plot([0 0],[ymin ymax],'k--')
plot([xmin 0],[0 0],'k--')
plot([xmin xmax],[1.5 1.5],'k--')
%plot([xmin xmax],[GMSTC2022 GMSTC2022],'k--')
hold off
flab = gca;
flab.XMinorTick = 'on';
flab.YMinorTick = 'on';
xtickformat('%.0f');
ytickformat('%.1f');
%flab.XTick = [xmin:xinc:xmax];
%flab.XTickLabel = [xmin:xinc:xmax];
%flab.YTick = [1 10 100 1000 10000];
%flab.YTickLabel = '';%[10 10 10 10 10];
flab.FontSize = 30;
flab.FontName = 'helvetica';
flab.FontWeight = 'normal';
%xlabel('accumulated emissions (GtC)','fontsize',40,'fontname','helvetica')
%ylabel('temperature change (    )','fontsize',40,'fontname','helvetica')
box('on')
axis([xmin xmax ymin ymax])
axis('square')

%-----------Forster bar chart ---------------------
%compare against Table S4 in https://essd.copernicus.org/articles/16/2625/2024/essd-16-2625-2024-supplement.pdf
M = [1.38 1.21 1.55  %ROF
1.29 1.16 1.40       %GWI
1.26 1.10 1.42       %KCC
1.31 1.10 1.60      %Forster et al., 2024 ensemble
HIW2023Bmu HIW2023Bmu-2*HIW2023Bsig HIW2023Bmu+2*HIW2023Bsig;
HIW2023mu HIW2023mu-2*HIW2023sig HIW2023mu+2*HIW2023sig];

xmin = 0;xmax = 7;xinc = 1;
ymin = 0;ymax = 2;yinc = 0.2;
figure(4);clf
hold on
f1 = bar((1:4)',M(1:4,1)); 
f2 = bar((5:6)',M(5:6,1)); 
n = 1;line([n n],M(n,2:3),'color',[0 0 0]);line([n-.1 n+.1],[M(n,2) M(n,2)],'color',[0 0 0]);line([n-.1 n+.1],[M(n,3) M(n,3)],'color',[0 0 0])
n = 2;line([n n],M(n,2:3),'color',[0 0 0]);line([n-.1 n+.1],[M(n,2) M(n,2)],'color',[0 0 0]);line([n-.1 n+.1],[M(n,3) M(n,3)],'color',[0 0 0])
n = 3;line([n n],M(n,2:3),'color',[0 0 0]);line([n-.1 n+.1],[M(n,2) M(n,2)],'color',[0 0 0]);line([n-.1 n+.1],[M(n,3) M(n,3)],'color',[0 0 0])
n = 4;line([n n],M(n,2:3),'color',[0 0 0]);line([n-.1 n+.1],[M(n,2) M(n,2)],'color',[0 0 0]);line([n-.1 n+.1],[M(n,3) M(n,3)],'color',[0 0 0])
n = 5;line([n n],M(n,2:3),'color',[0 0 0]);line([n-.1 n+.1],[M(n,2) M(n,2)],'color',[0 0 0]);line([n-.1 n+.1],[M(n,3) M(n,3)],'color',[0 0 0])
n = 6;line([n n],M(n,2:3),'color',[0 0 0]);line([n-.1 n+.1],[M(n,2) M(n,2)],'color',[0 0 0]);line([n-.1 n+.1],[M(n,3) M(n,3)],'color',[0 0 0])
plot([0 7],[1.5 1.5],'k--')
hold off
set(f1,'facecolor',[1 1 1]*.5,'edgecolor',[1 1 1]*.5,'facealpha',0.35,'edgealpha',0.35);
set(f2,'facecolor',COL,'edgecolor',COL,'facealpha',0.75,'edgealpha',0.75);
flab = gca; 
flab.XMinorTick = 'off';
flab.YMinorTick = 'on';
xtickformat('%.0f');
ytickformat('%.1f');
flab.XTick = '';
flab.XTickLabel = '';
flab.YTickLabel = '';
flab.FontSize = 30;
flab.FontName = 'helvetica';
flab.FontWeight = 'normal';
box('on')
axis([xmin xmax ymin ymax])
axis('square')

%------------------------------
%Fig.2. Test for persistence of 20 year out of sample forecast for
%sensitivity 
xmin = 1900-2;xmax = 2027;xinc = 10;
ymin = -2;ymax = 3;yinc = 0.2;
%recursive WLS estimates up to 2003, (y-c)/x thereafter
X = RECURSIVERESULT;
i = find(X(:,1)<=2003);
j = find(X(:,1)>=2003);
figure(5);clf
hold on
%plot recursive m
plot(X(i,1),X(i,2)*100,'col',COL*.8)
h = patch([X(i,1);flipud(X(i,1))],[(X(i,2)-2*X(i,3));flipud((X(i,2)+2*X(i,3)))]*100,[1 1 1]);
set(h,'facecolor',[1 0 0]*.5,'edgecolor',[1 0 0]*.5,'facealpha',0.25,'edgealpha',0.25);
%plot recursive c
plot(X(i,1),X(i,4)*2,'k')
h = patch([X(i,1);flipud(X(i,1))],[(X(i,4)-2*X(i,5));flipud((X(i,4)+2*X(i,5)))]*2,[1 1 1]);
set(h,'facecolor',[1 1 1]*.5,'edgecolor',[1 1 1]*.5,'facealpha',0.35,'edgealpha',0.35);
%2023 error bars
plot([2025 2025],[MM-2*MMsig MM+2*MMsig]*100,'col',COL*.8)
plot([2024 2026],[MM+2*MMsig MM+2*MMsig]*100,'col',COL*.8)
plot([2024 2026],[MM-2*MMsig MM-2*MMsig]*100,'col',COL*.8)
plot([2025 2025],[CC-2*CCsig CC+2*CCsig]*2,'k')
plot([2024 2026],[CC+2*CCsig CC+2*CCsig]*2,'k')
plot([2024 2026],[CC-2*CCsig CC-2*CCsig]*2,'k')
%observed sensitivity
i = find(yearMCS>=2003 & yearMCS<=2023);
h = patch([yearMCS(i);flipud(yearMCS(i))],[SO05(i);flipud(SO95(i))]*100,[1 1 1]);
set(h,'facecolor',[1 0 0]*.5,'edgecolor',[1 0 0]*.5,'facealpha',0.4,'edgealpha',0.4);
%zero lines
plot([xmin xmax],[0 0],'k:')
plot([2003 2003],[ymin ymax],'k:')
hold off
flab = gca;
flab.XMinorTick = 'off';
flab.YMinorTick = 'on';
xtickformat('%.0f');
ytickformat('%.1f');
flab.FontSize = 30;
flab.FontName = 'helvetica';
flab.FontWeight = 'normal';
box('on')
axis([xmin xmax ymin ymax])
axis('square')

return
%--------------------END-----------------------