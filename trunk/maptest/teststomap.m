function [mr,mr0,vlist] = teststomap(A)


n2 = size(A,1);
n  = fix(sqrt(n2));
[M,x2] = rdwalkM(n,[1,1,1]);
Iden = speye(n2);

nmax  = 500;
mr    = [];
mr0   = [];
vlist = [];

row1  = sparse(1,n^2);
row1(1) = 1;



for i = 0:10:nmax
   Ap    = @(x) probMap(A,M,i,x);
   [V,D] = eigs(Ap,n2);

   d     = sort(abs(diag(D)));
   mr    = [mr d];

   

   vlist = [vlist row1*x2];    
   row1  = row1*M;

   Ap    = @(x) probMap(Iden,M,i,x);
   [V0,D0] = eigs(Ap,n2);
   d0     = sort(abs(diag(D0)));
   mr0    = [mr0 d0]; 

end

plot(vlist,(1-mr)./(1-mr0));

