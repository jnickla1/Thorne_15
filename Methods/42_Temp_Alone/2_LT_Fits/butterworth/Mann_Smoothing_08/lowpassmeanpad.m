function [smoothed,mse] = lowpassmeanpad(indata,frequency)
%
% [smoothed,mse] = lowpassmeanpad(indata,frequency)
% lowpass filter.
% assumes unit time spaced series 
% 
% lowpass cutoff at f="frequency" (in cycles/time unit)
%
% The lowpass routine employs a 10 point butterworth filter. 
%
% pad either boundary w/ mean over last 1/2 filter width
% 
%
ipts=10;
fn=frequency*2;
npad=3*round(1/fn);
npad0=round(1/fn);
nn=length(indata);
padded(npad+1:npad+nn)=indata;
padded(npad+nn+1:nn+2*npad)=mean(indata(nn-npad0+1:nn));
padded(1:npad)=mean(indata(1:npad0));
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





