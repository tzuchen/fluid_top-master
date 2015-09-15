function sigma = mixnormI(yb,zb,n,method)
% Given the boundry between Red and Blue liquid,
% This function approximates the mix norm:
%   sigma =  <(I-<I>)^2>^(0.5)
% by area integral. 
% 
%
% Tzu-Chen Liang 2-7-2005

switch method

case 0
 p         = 10;
 nresol    = n*p;
 linespace = (0.5:nresol)/nresol;
 [ptX,ptY] = meshgrid(linespace,linespace);
 In = inpolygon(ptX,ptY,[0 yb 1],[0 zb 0]);
 Iave = sum(sum(In))/nresol^2;
 Itotal = 0;
 
 for i = 1:n
   for j = 1:n
      Iij = sum(sum(In((i-1)*p+1:i*p,(j-1)*p+1:j*p)))/p/p; 
      Iv  = (Iij-Iave)^2;
      Itotal = Itotal+Iv; 
   end
end

 sigma = (Itotal / n^2)^0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alternative approach
case 1

 lsp = (0:n)/n;
 %Iave   = dblquad(@(x,y) fun(x,y,yb,zb),0,1,0,1);
 Itotal = 0;
 I = zeros(n,n); 

 for i = 1:n
   for j = 1:n     
      I(i,j) = dblquad(@(x,y) fun(x,y,yb,zb),lsp(i),lsp(i+1),lsp(j),lsp(j+1))*n*n;
      %Iv  = (Iij-Iave)^2;
      %Itotal = Itotal+Iv; 
   end
 end
 
 %Iave
 Iave = sum(sum(I))/n/n;
 sigma = (sum(sum((I-Iave).^2))/n/n)^0.5;


 %sigma = (Itotal / n^2)^0.5;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function in = fun(x,y,yb,zb) 
 y = y*ones(length(x),1); 
 in = inpolygon(x,y,[0 yb 1],[0 zb 0]);














