% Copyright (c) 2016.
% All rights reserved. Please read the 'license.txt' for license terms.
% 
% Developers: Zhen Zhang, Pakorn Kanchanawong
% Contact: biekp@nus.edu.sg
function flag =  LargeCurFlag(x, y) % coordinates consistent with image in Matlab
x = x-min(x)+1;
y = y-min(y)+1;

temp = zeros(max(x)-min(x)+3, max(y)-min(y)+3);

x = x+1;
y = y+1;
temp(sub2ind(size(temp),x,y)) = 1;

s1 = temp(sub2ind(size(temp),x-1,y-1));
s2 = temp(sub2ind(size(temp),x-1,y));
s3 = temp(sub2ind(size(temp),x-1,y+1));
s4 = temp(sub2ind(size(temp),x,y-1));
s5 = temp(sub2ind(size(temp),x,y));
s6 = temp(sub2ind(size(temp),x,y+1));
s7 = temp(sub2ind(size(temp),x+1,y-1));
s8 = temp(sub2ind(size(temp),x+1,y));
s9 = temp(sub2ind(size(temp),x+1,y+1));

s = s1+s2+s3+s4+s5+s6+s7+s8+s9;
tipidx = find(s==2);

if length(tipidx)~=2
    flag = 1;
else
    x_change = x;
    y_change = y;
    tipchange = [x(tipidx(1))  y(tipidx(1))];
    d = 0;
    for i = 1:length(x)
        k = dsearchn([x_change  y_change],tipchange);
        d = d + pdist([tipchange; [x_change(k)  y_change(k)]]);
        
        tipchange = [x_change(k)  y_change(k)];
        x_change(k) = [];
        y_change(k) = [];
    end
    
    if 1.2 * pdist([[x(tipidx(1))  y(tipidx(1))];[x(tipidx(2))  y(tipidx(2))]])<d
        flag = 1;
    else
        flag = 0;
    end
    
end


