function mu = mun(n)

mu =1;
 for i =1:n-1
    mu = Si(mu);
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = Si(x)
y = (1-sqrt(1-x))/2;
