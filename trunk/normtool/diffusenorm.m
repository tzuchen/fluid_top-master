function [w] = diffusenorm(X,domain,d)
%
% This function evaluate the diffuse norm
% of X. with diffusion rate d.
%
%
% Tzu-Chen Liang 7-10-2006

if prod(size(X)) ==1
% return the weighting
  n = X;
  w = exp(-2*4*pi^2*d*(0:n).^2); 
else
   n = 2;
   k = size(X,1);
   if domain == 'f'
      Xf = X;
   elseif domain == 't'
      Xf = fft2(X);
   else
      error('domain must be t or f')
   end

   nlist =  [k/2:-1:1,0:1:k/2-1];

   P     =  ones(k,1)*nlist.^2+ (ones(k,1)*nlist.^2)';
   M     =  exp(-2*4*pi^2*d*P);
   M     =  fftshift(M);
   w     =  sqrt(sum(sum(M.*(abs(Xf).^2))))/k/k;
end
