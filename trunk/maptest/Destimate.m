function D = Destimate(mapfunction,param,n)
% Given a map and number of grids
% This function esitmates the effective D*
% 
% Tzu-Chen Liang 7-28-2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    wn     = n/10;
    dy     = 1/n;
    x      = dy/2:dy:1-dy/2;
    y      = dy/2:dy:1-dy/2;
    [X,Y]  = meshgrid(x,y);

  [Astruc] = maprefine2(n,[],mapfunction,param{:});

   [Xp,Yp] = feval(mapfunction,X,Y,1,1,param);
   Z       = sin(2*pi*Yp*wn)'.*sin(2*pi*Xp*wn)';
   z       = reshape(Z,n^2,1);
   z2      = Astruc.A*z;
   D       = -log(norm(z2)/norm(z))/(4*pi^2*2*wn^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
