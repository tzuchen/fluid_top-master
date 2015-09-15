function x = sdppotential(X,S,n,lo)


     x = (n+lo)*log(trace(X*S))-log(det(X)*det(S));


 if imag(x)~=0

   x = inf;
end
