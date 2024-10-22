function[smoothed0,icb,ice,mse0] = lowpassmin(indata,frequency)
%
% [smoothed,icb,ice,mse0] = lowpassmin(indata,frequency)
mse0 = 999.0;
for iconb=0:2
   for icone=0:2
       [smoothed,mse] = lowpass(indata,frequency,iconb,icone);
       if (mse<=mse0) 
           icb=iconb;
           ice=icone;
           mse0=mse;
           smoothed0=smoothed;
       end
   end
end





