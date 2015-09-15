function [Mstruc] = fourierdiffuseM(n,k,dt,period)
%
% This function generate a nxn matrix M
% To diffuse X, use the following command 
%
%  X   = real(ifft2(fft2(X).*M));
%
%  Currently, only periodical M..
%  k and dt are calibrated by rdwalkM 
%
%    Tzu-Chen Liang 6-27-2006


%ke      = k *1000;
%dte     = dt*1000; 

%dx      = 1/n;
%[Xs,Ys] = meshgrid(dx/2:dx:1-dx/2,dx/2:dx:1-dx/2);
%Xs      = Xs-0.5;
%Ys      = Ys-0.5;

%[th,r]  = cart2pol(Xs,Ys);
%df      = r.^2;

%M       = exp(-pi/2*ke*dte*df);
%M       = M/max(max(M));
%surf(M)
%M       = fftshift(M);


 nlist =  [n/2:-1:1,0:1:n/2-1];

 P     =  ones(n,1)*nlist.^2+ (ones(n,1)*nlist.^2)';
 M     =  exp(-4*pi^2*k*dt*P);
 M     =  fftshift(M);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Mstruc.M = M;
Mstruc.n      = n;
Mstruc.k      = k;
Mstruc.dt     = dt;
Mstruc.period = period;
