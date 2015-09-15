function [varyc,A,ycom] = ndcubesimu(n)

m = 2^n;

A = sparse(m,m);
ndig = length(dec2bin(m-1));


x= dec2bin(0:m-1,ndig);

for i = 1:m
  for j = 1:m
    if sum(abs(x(i,:) -x(j,:)))<=1;
      A(i,j)=1/(n+1);
    end
  end
end

[V,D] = eig(full(A'));



y = zeros(m,1);
y(1) = 1;
yf   = 1/m*ones(m,1);
yc = V'*y;
yfc  = V'*yf;

for i = 1:10
   vary(i)=sum(abs(y-yf)); 
   %varyc(i) =  sum(abs(yc-yfc));
   varyc(i) = sqrt(sum(((y./yf-1).^2).*yf))
   ycom(i,:) = (yc-yfc);

   yp = y;
   y = A'*y;
  % yc = D*yc;
 %  V'*(y-yf)-D*V'*(yp-yf)
 
 % norm(y-yf,1)-norm(yc-yfc,1)
  
 

% norm(y./yf-1,2)-norm(V'*(y./yf-1),2)
  

 


end


close all
figure
hold on
%plot(vary)
plot(varyc,'r')
%plot(sum(abs(ycom')),'g')

%yfc
%diag(D)
