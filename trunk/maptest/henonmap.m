function [xp,yp] = henonmap(x,y,period,i,param)
% functoion for henon map
% period = 1, periodical with period 1
%          0, not periodical
% i      > 0, henon map
%        < 0, inverse henon map
%
% param{1} = alpha
% param{2} = beta
% param{3} = modulus
%
%  Tzu-Chen Liang 3-24-2006
%



if nargin <5
   alpha = 0.2;
   beta  = -0.9999;
   modn  = 1;
else
   alpha = param{1};
   beta  = param{2};
   modn  = param{3};
end 

x = x * modn;
y = y * modn;

if i > 0
   if period ==1
      xp = mod(1-alpha*x.^2+y,modn);
      yp = mod(beta*x,modn);
   else
      xp =  1-alpha*x.^2+y;
      yp =  beta*x;
   end
else

   if period ==1
      xp = mod(y/beta,modn);
      yp = mod(alpha/beta^2*y.^2+x-1,modn);
   else
      xp = y/beta;
      yp = alpha/beta^2*y.^2+x-1;
   end
end

xp = xp / modn;
yp = yp / modn;

