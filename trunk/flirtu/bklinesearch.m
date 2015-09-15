function t = bklinesearch(X,S,dX,dS,alpha,beta,n,lo)


t  = 1;
while sdppotential(X+t*dX,S+t*dS,n,lo) > sdppotential(X,S,n,lo)+alpha*t*dsdppotential(X,S,n,lo,dX,dS)
       t = beta*t;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gX = dsdppotential(X,S,n,lo,dX,dS)

  p =(n+lo)/(trace(X*S));
  gX = trace((p*S-inv(X))*dX) +trace((p*X-inv(S))*dS);
