function [S4,alphay] = rotation_frequency_dependence_Lanczos(p1,p2);

alphay = [0:90]';

% frequency dependence: Lanczos filter
ny1 = 6*12+1;
ny2 = 20*12+1;

npos = 0;
for i = 1:numel(alphay)
    command = sprintf('%d',i);
    fprintf([repmat('\b',1,npos) '%s'],command);
    npos = length(command);
    
    alpha = alphay(i);
    ap1 = cosd(alpha); ap2 = sind(alpha);
    
    p31 = ap1*p1 - ap2*p2;
    p41 = ap1*p2 + ap2*p1;
    
    p31h = lanfilt(p31,2,1/ny1,[],ny1);
    p31m = lanfilt(p31,1,[],1/ny2,ny2);
    p31l = p31 - p31h - p31m;
    
    p41h = lanfilt(p41,2,1/ny1,[],ny1);
    p41m = lanfilt(p41,1,[],1/ny2,ny2);
    p41l = p41 - p41h - p41m;
    
%     s21 = nancorr(p31h,p41h);
%     s22 = nancorr(p31l,p41l);
%     s23 = nancorr(p31m,p41m);
%     S4(i) = s21^2+s22^2+s23^2;

    s21 = nancov(p31h,p41h);
    s22 = nancov(p31l,p41l);
    s23 = nancov(p31m,p41m);
    
    S4(i) = s21(1,2)^2+s22(1,2)^2+s23(1,2)^2;
end