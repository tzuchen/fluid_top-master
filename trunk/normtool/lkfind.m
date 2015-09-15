function lk = lkfind(n,klist)
%  
%  This function calculate the weighting of each wave number
%  to evaluate the mix-norm. 
%  Because the numerical integration is required, we do this 
%  only once, and save the result as a variable lk2000.
%
% Tzu-Chen Liang 7-7-2006



warning off all
for i = 1:length(klist)

  if klist(i) ==0
    lk(i)  = 1;
  else
    lk(i)  = lambdak(n,klist(i));
  end
end
warning on all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x  = f(n,k,s)

x  = 2^(n-2)*n^2*gamma(n/2)^2 *(besselj(n/2,s*pi*k).^2)./(s*pi*k).^n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function lk  = lambdak(n,k)
  
  fnk = @(s) f(n,k,s); 
  
  lk = quad(fnk,0,1);
 

