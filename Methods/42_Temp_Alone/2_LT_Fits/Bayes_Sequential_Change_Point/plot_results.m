function [temps,uncer] = plot_results(parameters, x, X, Y, Py, P, k)
% This function draws samples from the exact posterior distribution on the
% number of change points, their locations, and the parameters of the
% regression model, and then plots the results.
set(0,'DefaultAxesTitleFontWeight','normal');

[N, m] =size(X);            % N = total # data points (N), m = # regressors
d_min = parameters(1);      % Minimum distance between adjacent change points    
k_0 = parameters(2);        % Hyperparameter for the prior on the regression coefficients
v_0 = parameters(3); sig_0 = parameters(4); % Hyperparameter for scaled inverse chi-square prior on the variance
k_max = parameters(5);      % Maximum number of change points
num_samp = parameters(6);   % Number of sampled solutions
offset=parameters(7);
beta0= zeros(m,1);          % Mean of multivariate normal prior on regression coefficients
I=eye(m);                   % m x m identity matrix


ap=0.0455; %2 sigma
ac=0.3173; %1 simga

cvalue=[26./255, 44./255, 105./255];
cstdevmean=[62./255, 207./255, 117./255];
cerrorbar=[52./255, 235./255, 235./255];
pcolor=[255./255,20./255,147./255]; %'deeppink';


% Variables used in sampling procedure:
samp_holder=zeros(num_samp,k_max);  % Contains each of the num_samp change point solutions
chgpt_loc=zeros(1,N);               % The number of times a change point is identified at each data point
BETA = zeros(m,N);                  % Holds the regression coefficients
confpreds=zeros(N,2);

for i=1:num_samp
    %******** (i) Sample a Number of Change Points ***********************
    num_chgpts = pick_k1(k)-1;  % Since we allow for 0 changepoints, function returns the index of the 'k' vector, 
                                % which is offset from the number of change points by 1 
    if(num_chgpts>0)
        
        %******** (ii) Sample the Location of the Change Points ***********
        back=N;
        for kk=num_chgpts:-1:2      % Start at the end of the time series and work backwards
            temp=zeros(1,back-1);
            for v=1:back-1
                temp(v)= P(kk-1,v)+Py(v+1,back);  % Build the vector to sample from
            end
            % TO AVOID UNDERFLOW, USE:
            % M_temp = max(temp); temp = temp - M_temp;
            total=log(sum(exp(temp)));  
            temp(:)=exp(temp(:)-total);  % Normalize the vector
            changepoint=pick_k1(temp);   % Sample the location of the change point
            chgpt_loc(changepoint)= chgpt_loc(changepoint) +1; % Keep track of change point locations
            samp_holder(i,kk)=changepoint;  % Keep track of individual change point solutions
            
            %******** (iii) Sample the Regression Parameters ***********
            % Regression Coefficients (Beta)
            XTy=X(changepoint+1:back,:)'*Y(changepoint+1:back);
            J=k_0*I+X(changepoint+1:back,:)'*X(changepoint+1:back,:);
            beta_hat=J\(k_0*beta0+XTy);     %inv(J)*(k_0*beta0+XTy)
            

               BETA(:,changepoint+1:back) = BETA(:,changepoint+1:back) +repmat(beta_hat,1,back-changepoint)/num_samp;  
               confpreds(changepoint+1:back,:)=confpreds(changepoint+1:back,:)+ ...
                      pc_interval(X(changepoint+1:back,:),Y(changepoint+1:back),beta_hat,ac,ap)/num_samp;
            
            back=changepoint;   % Now work with the next segment
        end
        
        % The final changepoint
        %******** (ii) Sample the Location of the Change Points ***********
        kk=1;
        temp=zeros(1,back-1);
        for v=1:back-1
            temp(v)= Py(1,v)+Py(v+1,back); %Build the vector to sample from
        end
        % TO AVOID UNDERFLOW, USE:
        % M_temp = max(temp); temp = temp - M_temp;
        total=log(sum(exp(temp)));
        temp(:)=exp(temp(:)-total); % Normalize the vector
        changepoint=pick_k1(temp);  % Sample the location of the change point
        chgpt_loc(changepoint)= chgpt_loc(changepoint) +1; % Keep track of change point locations
        samp_holder(i,kk)=changepoint;  % Keep track of individual change point solutions
            
        %******** (iii) Sample the Regression Parameters ***********
        % Regression Coefficients (Beta)
        XTy=X(changepoint+1:back,:)'*Y(changepoint+1:back); 
        J=k_0*I+X(changepoint+1:back,:)'*X(changepoint+1:back,:);
        beta_hat=J\(k_0*beta0+XTy);     %inv(J)*(k_0*beta0+XTy)
        

               BETA(:,changepoint+1:back) = BETA(:,changepoint+1:back) +repmat(beta_hat,1,back-changepoint)/num_samp;  
               confpreds(changepoint+1:back,:)=confpreds(changepoint+1:back,:)+ ...
                      pc_interval(X(changepoint+1:back,:),Y(changepoint+1:back),beta_hat,ac,ap)/num_samp;
        
        %The final sub-interval
        XTy=X(1:changepoint,:)'*Y(1:changepoint); 
        J=k_0*I+X(1:changepoint,:)'*X(1:changepoint,:); 
        beta_hat=J\(k_0*beta0+XTy);     %inv(J)*(k_0*beta0+XTy)
        

               BETA(:,1:changepoint) = BETA(:,1:changepoint) +repmat(beta_hat,1,changepoint)/num_samp;  
               confpreds(1:changepoint,:)=confpreds(1:changepoint,:)+ ...
                      pc_interval(X(1:changepoint,:),Y(1:changepoint),beta_hat,ac,ap)/num_samp;

    else    % 0 change points, so a single homogeneous segment
        XTy=X'*Y; 
        J=k_0*I+X'*X; 
        
        %******** (iii) Sample the Regression Parameters ***********
        % Regression Coefficients (Beta)
        beta_hat=J\(k_0*beta0+XTy);     %inv(J)*(k_0*beta0+XTy)

               BETA = BETA +repmat(beta_hat,1,N)/num_samp;
               confpreds=confpreds+ ...
                      pc_interval(X,Y,beta_hat,ac,ap)/num_samp;

    end
end

%BETA=BETA/num_samp;             % Average regression coefficient at each data point.
chgpt_loc=chgpt_loc/num_samp;   % Posterior probability of a change point

%**********(iv) Plot the Results ********************************
% Adapt as Necessary - Axes are specific to temperature anomalies data set

%start
%end
%
%n=end-start+1;

post_lim=0.8;

model=zeros(1,N);
for i=1:N
    model(i)=X(i,:)*BETA(:,i)+offset;
end


%figure(1); 
%clf; 
%hold

%[ax, h1, h2] = plotyy(x,model,x,chgpt_loc, 'plot');

%inBetweenstd = [model'-confpreds(:,2);  ...
%    flipud(model'+confpreds(:,2))];
%fill([x; flipud(x)], inBetweenstd, cerrorbar , 'LineStyle','none');
%inBetweenstm = [model'-confpreds(:,1)*2;  ...
%    flipud(model'+confpreds(:,1)*2)];
%fill([x; flipud(x)], inBetweenstm, cstdevmean , 'LineStyle','none');

%plot(x,Y+offset,'.','MarkerSize',8, 'MarkerEdgeColor',[0.5,0.5,0.5]);



%title('Change Point Lines, sampled ' +string(num_samp)+' times', 'fontsize', 12)

%M=cat(1,x',model,confpreds(:,1)',confpreds(:,2)',chgpt_loc);
%csvwrite('BSCtemperatures.csv',M);
temps = model;
uncer = confpreds(:,1)';

%gca=ax(1);
%set(gca, 'Xlim', [1850 2025])
%set(gca, 'Xtick', 1850:25:2025)
%set(gca, 'Xticklabel', 1850:25:2025)
%xlabel('Year')
%set(gca, 'Ylim', [-1.0 2.5])
%set(gca, 'Ytick', -1.0:0.2:2.5)

%set(gca, 'Yticklabel', -1.0:0.2:2.5)
%set(gca.YAxis, 'TickLabelFormat','%.1f')
%ytickformat('%.1f')

%set(gca,'Ycolor', 'k')      % Left Y-axis colored black
%set(gca, 'fontsize', 12)    % Default is 10
%set(get(ax(1), 'Ylabel'), 'String', 'Temperature (K)', 'fontsize', 12)
%set(h1, 'Color', cvalue)       % Model plotted in green, default blue
%set(h1, 'LineStyle', '-'); % Default is a solid line
%set(h1, 'LineWidth', 2)     % Default is 0.5

%gca=ax(2);
%set(gca, 'Xlim', [1850 2025])
%set(gca, 'Xtick', 1850:25:2025)
%set(gca, 'Xticklabel', 1850:25:2025)
%set(gca, 'Ylim', [0 2.5*post_lim/10*11])
%ax(2).TickLength = [0.005,0.005];
%set(gca, 'Ycolor', pcolor)     % Right Y-Axis colored red
%set(gca, 'YDir' , 'reverse')
%set(gca, 'Yticklabel', 0:0.1:post_lim)
%set(gca, 'Ytick', 0:0.1:post_lim)
%set(gca, 'fontsize', 12)
%set(get(ax(2), 'Ylabel'), 'String', 'Posterior Probability', 'fontsize', 12)
%ylh = get(ax(2), 'Ylabel')
%ylh.Position(2) = ylh.Position(2) - abs(ylh.Position(2) * 0.65);
%set(h2, 'Color', pcolor)       % Posterior probabilities plotted in red
%set(h2, 'LineWidth', 2)
%line(ax(2), [1850 2025],[post_lim post_lim], 'Color',pcolor, 'LineWidth', 0.2);    
        %Draw a line representing the end of the data set

%line( [x(N); x(N)], [-0.8; 0.8], 'Color','k', 'LineWidth', 2);    
        %Draw a line representing the end of the data set

%hold;

%inset=axes('Position',[0.14,0.628,0.16,0.298]);
%box on
%bar(k,'LineStyle','none','FaceColor',pcolor)
%set(inset, 'Ycolor', pcolor)     % Right Y-Axis colored red
%set(inset, 'YDir' , 'reverse')
%set(inset, 'Ylim', [0 post_lim])
%set(inset, 'Yticklabel', 0:0.1:post_lim)
%set(inset, 'Ytick', 0:0.1:post_lim)
%set(inset, 'fontsize', 12)
%set(inset,'YAxisLocation','right')
%set(inset, 'YDir' , 'reverse')
%set(inset,'YAxisLocation','right')
%set(inset, 'Xlim', [0 8])
%xlabel('# of Lines')

%bring the post dist to front
%uistack(h1,'top')


end

function confpred=pc_interval(X,y,beta_hat,alpha_c, alpha_p)
%note m must be number of columns in X
    n=length(y);%
    m=length(beta_hat);%2
    model = X*beta_hat;
    %SE regresion line
    serl=sqrt(sum((y-model).^2)/(n-m));
    x_=mean(X(:,2));
    xsqdist=((X(:,2)-x_).^2);
    
    conf= -tinv(alpha_c/2,n-m)* serl* sqrt(1/n + xsqdist/sum(xsqdist));
    pred= -tinv(alpha_p/2,n-m)* serl* sqrt(1+ (1/n) + xsqdist/sum(xsqdist));
    confpred = [conf,pred];
end
