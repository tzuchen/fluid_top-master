function [M] = rdwalkM(n,p,k,t,period)
%
%  This function generates a (n^2)*(n^2) matrix
%    
%  period = 1  => periodical  
%         = 0  => not periodical
%  
%
%  Tzu-Chen Liang 6-9-2006

xl = 1;

dx  = xl/n;

x = dx/2:dx:xl-dx/2;
y = dx/2:dx:xl-dx/2;

[X,Y] = meshgrid(x,y);
Xv    = reshape(X,n^2,1);
Yv    = reshape(Y,n^2,1);

b = (n^2)/(2*pi^2)*(1 - exp(-4*pi^2*k*t/p));

probd = [b 1-2*b b];  

 probd = probd/sum(probd)/2;

if period == 1
  D = diag(sparse(probd(1)*ones(n-1,1)),1) +...
      diag(sparse(probd(3)*ones(n-1,1)),-1)+...
      diag(sparse(probd(2)*ones(n  ,1)));
  D(1,n) = probd(3);
  D(n,1) = probd(1);
else
  D = diag(sparse(probd(1)*ones(n-1,1)),1) +...
      diag(sparse(probd(3)*ones(n-1,1)),-1)+...
      diag(sparse(probd(2)*ones(n  ,1)));
  D(1,1) = D(1,1)+ probd(3);
  D(n,n) = D(n,n)+ probd(1);
end

 I      = speye(n);

 %X2  = zeros(n,n);
 %for i= 1: n
 %  for j = 1:n
 %     im = min(i,n+2-i);
 %     jm = min(j,n+2-j);

 %    X2(i,j) = 1/n/n*((im-1)^2+(jm-1)^2);
 %  end
 %end 
 %x2 = reshape(X2,n^2,1);

 M = kron(D, I) + kron(I,D);
 M = M^p;

