function [S1,alphay] = rotation_spatial_variance(sp1,sp2,lons,lats);

alphay = [0:90]';
S1 = nan(numel(alphay),1);

% spatial variance
npos = 0;
for i = 1:numel(alphay)
    command = sprintf('%d',i);
    fprintf([repmat('\b',1,npos) '%s'],command);
    npos = length(command);
    
    alpha = alphay(i);
    ap1 = cosd(alpha); ap2 = sind(alpha);
   
    sp21 = ap1*sp1 - ap2*sp2;
    sp22 = ap1*sp2 + ap2*sp1;
    
    sp21_2 = sp21.^2;
    sp22_2 = sp22.^2;
    
    s = sp21_2.*sp22_2;
    S1(i) = area_weighted_mean(s,lons,lats);
end
