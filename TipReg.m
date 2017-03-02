% Copyright (c) 2016.
% All rights reserved. Please read the 'license.txt' for license terms.
% 
% Developers: Zhen Zhang, Pakorn Kanchanawong
% Contact: biekp@nus.edu.sg

function RegTip = TipReg(LL, label, curr_tip, RegR)
[x y] = find(LL==label);
xx = x-min(x)+1;
yy = y-min(y)+1;
tip = curr_tip(1,1:2)-[min(x) min(y)]+[1 1];

Dlist = sqrt((xx-tip(1)).^2 + (yy-tip(2)).^2);
xx1 = xx(find(Dlist<=RegR));
yy1 = yy(find(Dlist<=RegR));
temp = zeros(max(xx),max(yy));
temp(sub2ind(size(temp),xx1,yy1))=1;
center = regionprops(temp,'centroid');
center  = center(1).Centroid; % this has been confirmed that the result is [ Colume# (from left to right)   Row# (from top to bottom) ]
center  = [center(2)  center(1)];
RegTip(1,1) = check_dir_between2pts(center,tip);
RegTip(1,2:3) = center + [min(x)  min(y)];
end

