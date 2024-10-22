function AveragingComparison
close all
set(0,'DefaultAxesTitleFontWeight','normal');
data = readtable('../HadCRUT5.global.annual.csv');
temps=data.Var2+287;
dates=data.Var1;

cstdevmean=[62./255, 207./255, 117./255];
cerrorbar=[52./255, 235./255, 235./255];
cvalue=[26./255, 44./255, 105./255];
grey=[0.5,0.5,0.5];

N = 30;
moving_aves=zeros(length(temps),1);
std_aves=zeros(length(temps),1);
cN=ceil(N/2);
fN=floor(N/2);
fg = figure(1);
for i = (cN):(length(temps)-fN)
        lasta=i+fN;
        firsta=i-cN;
        moving_aves(i) = mean(temps((firsta+1):lasta));
        std_aves(i)= std(temps((firsta+1):lasta));
end


subplot(2,2,1)
setup_plot()
updN=10;
title(string(N)+'-year Mean, Updated Every '+string(updN)+' years')
for yea = (cN):updN:(length(temps)-fN)
    year=yea+1850;
    rectangle('Position', ...
        [(year-cN+1) (moving_aves(yea)-std_aves(yea)/sqrt(N)) N std_aves(yea)/sqrt(N)*2], ...
        'LineStyle', 'none', 'FaceColor', cstdevmean);
    line([(year-cN+1) (year+fN+0.9)],[moving_aves(yea) moving_aves(yea)],'LineWidth',1, 'Color',cvalue);
    e = errorbar(year,moving_aves(yea),std_aves(yea)*2, 'LineWidth',2,'Marker','none', 'Color',cerrorbar );
end
plot(dates,temps, '.', 'MarkerEdgeColor',grey);


subplot(2,2,2)
hold on
rax=dates(cN:(length(temps)-fN));
ray=moving_aves(cN:(length(temps)-fN));
inBetweenstd = [ray-std_aves(cN:(length(temps)-fN))*2;  ...
    flipud(ray+std_aves(cN:(length(temps)-fN))*2)];
fill([rax; flipud(rax)], inBetweenstd, cerrorbar , 'LineStyle','none');
inBetweenstm = [ray-std_aves(cN:(length(temps)-fN))/sqrt(N)*2;  ...
    flipud(ray+std_aves(cN:(length(temps)-fN))/sqrt(N)*2)];
fill([rax; flipud(rax)], inBetweenstm, cstdevmean , 'LineStyle','none');
plot(rax,ray, 'Marker', 'none','LineWidth',1,'Color', cvalue);
plot(dates,temps, '.', 'MarkerEdgeColor',grey);
title('Running '+string(N)+'-year Mean')
setup_plot()
plot(dates,temps, '.', 'MarkerEdgeColor',grey);



subplot(2,2,3)
setup_plot()
title('Optimal Climate Normal Chunks')
et50=1850-1;
computing=2;
Nguess1=30;
sp=1; %start
ep=Nguess1; %end
tau=0;
Noptimal=30;
Nguess=Nguess1;
OCNM=[];

while computing>=1
    
    %generate trend, sigma, beta, g for this guess interval Nguess= 30
    p = polyfit(dates(sp:ep),temps(sp:ep),1) ;
    model = p(1)*dates(sp:ep)+ p(2)*ones(ep-sp+1,1);
    sigma = sqrt(sum((model-temps(sp:ep)).^2)/(length(model)-2)); % regression std error
    beta = p(1)/sigma;
    g = (-model(1:end-1)+temps(sp:ep-1))'* (-model(2:end)+temps(sp+1:ep))  /sum((model-temps(sp:ep)).^2);
    plot(dates(sp:ep),model, 'Marker', 'none','LineWidth',0.5,'Color', grey);
    
    %find Noptimal
    funct = @(N)(OCNerr(N,g,beta,tau+.1));
    Noptimal = fminbnd(funct,1,50);
    if Noptimal<2
        Noptimal=2;
    end
    if round(Noptimal)+sp>length(temps)
        Noptimal=length(temps)-sp;
        computing=0;% extended too far
    end
    if Noptimal<=Nguess1
        Nguess=round((Noptimal+Nguess1)/2);
    else
        Nguess=Noptimal;
        continue
    end
        
    %compute stats and plot for start_p through start_p+Noptimal-1
    ep=sp+round(Noptimal)-1;
    if ep<=length(temps)
        avg_now=mean(temps(sp:ep));
        std_now=std(temps(sp:ep));
    else
        real_data=length(temps)-sp;
        proj_data=ep-length(temps);
        temps2=[temps(sp:end);p(1)*[dates(end)+1:dates(end)+proj_data]'+p(2)*ones(1,proj_data)];
        avg_now=mean(temps2(sp:ep));
        std_now=std(temps2(sp:ep));
    end
    
    rectangle('Position', ...
        [sp+et50 (avg_now-std_now/sqrt(Noptimal)) Noptimal std_now/sqrt(Noptimal)*2], ...
        'LineStyle', 'none', 'FaceColor', cstdevmean);
    line([(sp+et50) (ep+1+et50)],[avg_now avg_now],'LineWidth',1, 'Color',cvalue);
    e = errorbar(sp+Noptimal/2+et50,avg_now,std_now*2, 'LineWidth',2,'Marker','none', 'Color',cerrorbar );
    OCNM=[OCNM,[[(sp+et50):(ep+et50)];repmat([avg_now; std_now/sqrt(Noptimal);std_now*2],1,round(Noptimal))]];
    
    
    %reset interval and such
    sp=sp+round(Noptimal);
    if sp+1>=length(temps)
        computing=0;
        break
    end

    
    if sp+Nguess>length(temps)
        ep=length(temps);
        computing=1;
        %tau=Nguess-ep+sp;
    else
        ep=sp+Nguess;
    end
        
end

csvwrite('OCNtemperatures.csv',OCNM);

plot(dates,temps, '.', 'MarkerEdgeColor',grey);


subplot(2,2,4)
setup_plot()
title('Butterworth Smoothed')
smoothederrorbar=lowpassmeanpad(temps,1/3);
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

AddLetters2Plots({'a)', 'b)', 'c)','d)'})

buswrite=[dates,smoothedavg,maxs-mins,smoothedavg-smoothederrorbar];
csvwrite('BUStemperatures.csv',buswrite');

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

function nerr=OCNerr(N,g,b,t)
    nerr= (1+g)/(1+g+(N-1)*(1-g));
    nerr = nerr + (b*((N-1)/2+t))^2;
end
