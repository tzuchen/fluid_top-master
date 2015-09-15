function [xf,yf] = doubletwistmap(x,y,period,i,param)
% 
% This function calls @twistmap and creates
% a double twist map. 
%  
% 
% Tzu-Chen Liang   6-20-2006

paramat = cell2mat(param);
param1   = num2cell(paramat(1:3));
param2   = num2cell(paramat(4:6));

xc1 = param1{1};
yc1 = param1{2};
xc2 = param2{1};
yc2 = param2{2};
 

[xf,yf] = twistmap(x,y,period,i,param1);
dx1     = xf - x;
dy1     = yf - y;
[xf,yf] = twistmap(x,y,period,i,param2);
dx2     = xf - x;
dy2     = yf - y;


r1      = ((x-xc1).^2+(y-yc1).^2).^2;
r2      = ((x-xc2).^2+(y-yc2).^2).^2;
t       = r2./(r1+r2);
xf      = x + t.*dx1 + (1-t).*dx2;
yf      = y + t.*dy1 + (1-t).*dy2; 

