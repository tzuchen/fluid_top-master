
close all
clear k

n        = 1500;
param{1} = 0.2;
[Astruc] = maprefine2(n,[],@standardmap,param{1});
wnlist   = [1:150];

  dy     = 1/n;
  x      = dy/2:dy:1-dy/2;
  y      = dy/2:dy:1-dy/2;
  [X,Y]  = meshgrid(x,y);

for i = 1: length(wnlist);

   wn      = wnlist(i);
   [Xp,Yp] = standardmap(X,Y,1,1,param);
   Z       = sin(2*pi*Yp*wn)'.*sin(2*pi*Xp*wn)';

   [Xa1f,var2,var1,mnorm,w] = testAf(Astruc,[],Z,2,0,0,0);

%imagesc(Z)
%figure
%imagesc(Xa1f)

 k(i) = -log(var2(2)/var2(1))/(4*pi^2*2*wn^2);
end
 

  plot(wnlist,k)
