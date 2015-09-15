function [A,B] = maprefine1d3(n,mapfunction1d)
% 1-d version maprefine function
% the final correct version 
% Using the definition of koopman operator
%
% Tzu-Chen Liang   6-7-2007
%


dx=1/n;
x = 0:dx:1;
x(1) = 1e-8;
x(end) = 1-1e-9;

y = feval(mapfunction1d,x);

yg = ceil(y*n);
nrow = abs(diff(yg))+1;
nnzsize = sum(nrow);
iind = zeros(nnzsize,1);
jind = zeros(nnzsize,1);
s    = zeros(nnzsize,1);

ygmin = min(yg(1:n),yg(2:n+1));
ygmax = max(yg(1:n),yg(2:n+1));
yimin = min(y(1:n),y(2:n+1));
yimax = max(y(1:n),y(2:n+1));

p = 1;
count = 0;
for i = 1:n 


  iind(p:p+nrow(i)-1) = ygmin(i):ygmax(i);
  jind(p:p+nrow(i)-1) = i;
  
  s(p) = abs(yimin(i)/dx-(ygmin(i)));
  s(p+1:p+nrow(i)-2) = 1;
  s(p+nrow(i)-1) = abs(yimax(i)/dx-(ygmax(i)-1));


  p = p+nrow(i);
end




B = sparse(jind,iind,s,n,n);

A = sparse(diag(sparse(1./sum(B'))))*B;
B = B*sparse(diag(sparse(1./sum(B))));
B = B';


