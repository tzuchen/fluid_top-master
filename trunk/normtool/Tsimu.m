function nlist = Tsimu(T)


n  = size(T,1);
x0 = zeros(n,1);
x0(2)=1;
x  = x0;
nlist = norm(x);

for i = 1: 15
 x = T'*x;

 nlist=[nlist ,norm(x)]; 
end

nlist(end)/nlist(end-1)
semilogy(nlist);
