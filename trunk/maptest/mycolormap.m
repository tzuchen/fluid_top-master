function X = colormap(n,iter,mapfunction,param)
% This function uses the inverse map to simulate
% the color map.
% 
%  n = 500 is quite good
%
% Tzu-Chen Liang   3-24-2006
%
close all
figure 
dx      = 1/n;
xs      = dx/2 :dx:1-dx/2;
ys      = xs;
[Xs,Ys] = meshgrid(xs,ys);
xs      = reshape(Xs,n^2,1);
ys      = reshape(Ys,n^2,1);
par{1}  = param;

for i = 1:iter
    
    xva     = (cos(2*pi*xs));
    X       = reshape(xva,n,n);
    imagesc(X);
    M(i)    = getframe;
    [xs,ys] = feval(mapfunction,xs,ys,1,-1,par);
end










