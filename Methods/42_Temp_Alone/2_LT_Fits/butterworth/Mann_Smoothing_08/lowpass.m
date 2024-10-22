function [smoothed,mse] = lowpass(indata,frequency, iconb, icone)
%"minslp"
% [smoothed,mse] = lowpass(indata,frequency,iconb,icone)
% lowpass filter.
% assumes unit time spaced series 
% 
% lowpass cutoff at f="frequency" (in cycles/time unit)
%
% icon is any one of several possible boundary constraints on the smooth at 
% the begin/end x limits
% (0) minimum norm,
% (1) minimum slope
% (2) minimum roughness 
%
% The lowpass routine employs a 10 point butterworth filter. 

% Rather than implementing contraints (0)-(2) in the frequency domain (as in Park, 1992, 
% Ghil et al, 2003) we use the following approximate implementations of the boundary 
% constraints
% 
% (0) pad series with long-term mean value beyond x boundary
% (1) pad series with values over last 1/2 filter width reflected w.r.t. x
%     [imposes a local maximum/minimum at x boundary]
% (2) pad series with values over last 1/2 filter width reflected w.r.t. x and y (w.r.t. final value)
%     [imposes a point of inflection at x boundary]
%
ipts=10;
fn=frequency*2;
% pad far enough to remove any influence of boundaries of padded series on the interior
npad=3*round(1/fn);
nn=length(indata);
padded(npad+1:npad+nn)=indata;
padded(npad+nn+1:nn+2*npad)=indata(nn);
padded(1:npad)=indata(1);
for j=nn+npad+1:nn+2*npad
   ipad=j-nn-npad;
   if (icone==0) 
       apad=mean(indata);
   else if (icone==1)
           apad=indata(nn-ipad);
        else
           apad=2*indata(nn)-indata(nn-ipad);
        end
   end
   padded(j:j)=apad;
end
for j=1:npad
   ipad=j;
   if (iconb==0) 
       apad=mean(indata);
   else if (iconb==1)
           apad=indata(npad-ipad+1);
        else
           apad=2*indata(1)-indata(npad-ipad+1);
        end
   end
   padded(j:j)=apad;
end
%
% smoothed=padded;
[b,a]=butter(ipts,fn,'low');
smoothed0=filtfilt(b,a,padded);
smoothed=smoothed0(npad+1:nn+npad)';
%
% determine mse relative to raw data
%
asum = 0.0;
resid=smoothed-indata;
mse=var(resid)/var(indata);
