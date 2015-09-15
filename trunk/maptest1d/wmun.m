function w = wmun(y,k,n) 
%  Calculate w_(mu_n)^k
% 
%  Tzu-Chen Liang 5-9-2007



 mun =1;
 for i =1:n-1
  mun = Si(mun);
 end


 w  = wfun(y,k,mun);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = S(x)
y = 4*x.*(1-x);
function y = Si(x)
y = (1-sqrt(1-x))/2;
