function p = permutegenerate(n)

% np = n^3;
% nu = n^2*(n+1);
% 
% ind.u  = 1:nu;
% ind.v  = nu+1: 2*nu;
% ind.w  = 2*nu+1:3*nu;
% ind.p   = 3*nu+1:3*nu+np;
% 
% 
% P = reshape(ind.p,n,n^2);
% P = [ P ; zeros(1,n^2);];
% ind.p = reshape(P,1,nu);
% 
% p =reshape([ind.u; ind.v; ind.w; ind.p],4*nu,1);
% p = p(find(p>0));


np = n^3;
nu = n^2*(n-1);

ind.u  = 1:np;
ind.v  = np+1: np+nu;
ind.w  = nu+np+1:nu+nu+np;
ind.p  = nu+nu+np+1:nu+nu+np+np;


V = reshape(ind.v,n-1,n^2);
V = [ V ; zeros(1,n^2);];
ind.v = reshape(V,1,np);

W = reshape(ind.w,n-1,n^2);
W = [ W ; zeros(1,n^2);];
ind.w = reshape(W,1,np);


p =reshape([ind.u; ind.v; ind.w; ind.p],4*np,1);
p = p(find(p>0));
p = p(1:end-1);