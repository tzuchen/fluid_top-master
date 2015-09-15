function []=testinvariant(mdot,mapfunction)
% 
%  test whether mdot is an imvariant distribution of the mapfunction
%
%  Tzu-Chen Linag 3-24-2006


n       = fix(sqrt(length(mdot)));
Vm      = min(mdot);
VM      = max(mdot);
V       = reshape(mdot,n,n);
imagesc(V,[Vm,VM])

dx      = 1/n;
xs      = dx/2 :dx:1-dx/2;
ys      = xs;
[Xs,Ys] = meshgrid(xs,ys);
xs      = reshape(Xs,n^2,1);
ys      = reshape(Ys,n^2,1);

subplot(221)
scatter(xs,ys,30,mdot/Vm,'filled');
[xs,ys] = feval(mapfunction,xs,ys,1,-1);
subplot(222)
scatter(xs,ys,30,mdot/Vm,'filled');
[xs,ys] = feval(mapfunction,xs,ys,1,-1);
subplot(223)
scatter(xs,ys,30,mdot/Vm,'filled');
[xs,ys] = feval(mapfunction,xs,ys,1,-1);
subplot(224)
scatter(xs,ys,30,mdot/Vm,'filled');



