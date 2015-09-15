function A = maprefine1d2(n,mapfunction1d)
% 1-d version maprefine function
% the correct version 
%
% Tzu-Chen Liang   4-23-2007


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


p = 1;
count = 0;
for i = 1:n 
  iind(p:p+nrow(i)-1) = min(yg(i),yg(i+1)):max(yg(i),yg(i+1));
  jind(p:p+nrow(i)-1) = i;
  s(p:p+nrow(i)-1) = 1/nrow(i);
  p = p+nrow(i);
end




A = sparse(jind,iind,s,n,n);

%yd = diff(y);


