function w = mixnorm(X,domain,lk)
%
% This function calculate the mix-norm of X in R^nxn
% lk is the weighting of each wave number, it is found
% by the function lkfind(n,klist). 
%
%
% Tzu-Chen Liang 7-7-2006


 n = 2;
 k = size(X,1);
 if domain == 'f'
   Xf = X;
 elseif domain == 't'
   Xf = fft2(X);
 end

 
  nlist =  [k/2:-1:1,0:1:k/2-1];
  P     = (ones(k,1)*nlist.^2+ (ones(k,1)*nlist.^2)').^0.5;

  D = interp1(0:length(lk)-1,lk,P);
  D = fftshift(D);
  w = ((sum(sum(D.*abs(Xf).^2)))^0.5)/k/k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





