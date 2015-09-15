function X = createcoord(l,r,n,p)
% The function creates basic information of coordinates in 1 dimension.
% l: the left corner
% r: the right corner
% n: number of grids
% p: if it is periodical, p=1, otherwise p=0; 
% gisze: the grid size
%
% Tzu-Chen Liang   10/30/2005


X.l 	= l;
X.r 	= r;
X.n 	= n;
X.gsize = (r-l)/n;
X.p     = p;
