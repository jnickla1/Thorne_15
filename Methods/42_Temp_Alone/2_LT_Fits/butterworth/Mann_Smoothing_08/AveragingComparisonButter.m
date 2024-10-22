function AveragingComparisonButter
close all

set(0,'DefaultAxesTitleFontWeight','normal');
data = readtable('../HadCRUT5.global.annual.csv');
temps=data.Var2+287;
dates=data.Var1;
N = 30;

cstdevmean=[62./255, 207./255, 117./255];
cerrorbar=[52./255, 235./255, 235./255];
cvalue=[26./255, 44./255, 105./255];
grey=[0.5,0.5,0.5];

subplot(1,1,1)
setup_plot()
title('Butterworth Smoothed')
%for shift = 0:
smoothederrorbar=lowpassmeanpad(temps,1/5);
plot(dates,smoothederrorbar,'LineWidth',1, 'Color',cerrorbar);
smoothedavg=lowpassadaptive(temps,1/30);
mins=smoothedavg;
maxs=smoothedavg;
lout=3;
for yea = (50):N/2:(length(temps))
    ABC=[lowpassadaptive(temps(1:yea),1/30),lowpassmeanpad(temps(1:yea),1/30),lowpass(temps(1:yea),1/30,1,1)];
    plot(dates(1:yea),ABC(:,1),'LineWidth',1, 'Color',cstdevmean);
    plot(dates(1:yea),ABC(:,2),'LineWidth',1, 'Color',cstdevmean);
    plot(dates(1:yea),ABC(:,3),'LineWidth',1, 'Color',cstdevmean);
    if yea==length(temps)
        lout=0;
    end
    maxtemp = max([ABC,maxs(1:yea)],[],2);
    mintemp = min([ABC,mins(1:yea)],[],2);
    maxs(1:yea-lout)=mintemp(1:yea-lout);
    maxs(1:yea-lout)=maxtemp(1:yea-lout);
end
plot(dates,smoothedavg,'LineWidth',1, 'Color',cvalue);
plot(dates,temps, '.', 'MarkerEdgeColor',grey);
buswrite=[dates,smoothedavg,maxs-mins,smoothedavg-smoothederrorbar];

%csvwrite('BUStemperatures.csv',buswrite');
end

function setup_plot()
ylim([286.1 288.3])
xlim([1850,2025])
xticks(1850:25:2025)
yticks(286.2:0.2:288.2)
ytickformat('%.1f')
xlabel('Year')
ylabel('Temperature (K)')
box on
hold on
end
