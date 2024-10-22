function[smoothedbest,w0,w1,w2,msebest] = lowpassadaptive(indata,frequency)
%
% [smoothedbest,w0,w1,w2,msebest] = lowpassadaptive(indata,frequency)
%
% consider all possible weighted combinations of the three constraints
% (min norm/min slope/min roughness as approximately implemented 
% in routine 'lowpass)
% choose that which minimizes MSE
%
msebest = 1e+15;
% first optimize fixed lower constraint 
[smoothedlower,icb,ice,mse0] = lowpassmin(indata,frequency);
% no perform adaptive optimization of upper constraint
[smoothed0,mse0] = lowpass(indata,frequency,icb,0);
[smoothed1,mse1] = lowpass(indata,frequency,icb,1);
[smoothed2,mse2] = lowpass(indata,frequency,icb,2);
for weight0=0:0.01:1
 for weight1=0:0.01:1.0-weight0
       weight2=1-weight0-weight1;
       smoothed=weight0*smoothed0+weight1*smoothed1+weight2*smoothed2;
       mse = var(smoothed-indata)/var(indata);
       if (mse<=msebest) 
           w0=weight0;
           w1=weight1;
           w2=weight2;
           msebest=mse;
           smoothedbest=smoothed;
       end

 end
end





