function A = maprefine1d(n,mapfunction1d)
% 1-d version maprefine function
%
% Tzu-Chen Liang   4-18-2007


dx=1/n;
xlist = dx/2:dx:1-dx/2;

ylist = feval(mapfunction1d,xlist);

yn = ylist*n;

yl = yn - floor(yn);
yr = ceil(yn) - yn;


iind = [1:n, 1:n];
jind = [floor(yn), floor(yn)+1];
s    = [yl,yr];

A = sparse(iind,jind,s,n,n);
