function [xf,yf] = linkedtwistmap(x,y,period,i,param)
%
%  param{1} = center x
%  param{2} = center y
%  param{3} = rc  
%
%
% Linked Twisted Map
%
% Tzu-Chen Liang 6-13-2006

param1{1} = param{1};
param1{2} = param{2};
param1{3} = param{3};
param1{4} = param{4};

param2{1} = param{5};
param2{2} = param{6};
param2{3} = param{7};
param2{4} = param{8};


xf = x;
yf = y;

for i = 1:5
[xf,yf] = twistmap2(xf,yf,period,i,param1);
end
for i = 1:5
[xf,yf] = twistmap2(xf,yf,period,i,param2);
end
