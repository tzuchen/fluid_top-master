function N = Nestimate(mapfunction,param,D)
% Given a map and the D*
% This function esitmates number of grid n required 
% using linear interpolation
% 
% Tzu-Chen Liang 7-28-2006
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nl = 100;
nr = 500;

Dl = Destimate(mapfunction,param,nl);
Dr = Destimate(mapfunction,param,nr);



al =1/nl^2;
ar =1/nr^2;
ae = exp(log(al)+(log(D)-log(Dl))*(log(ar)-log(al))/(log(Dr)-log(Dl)));



N = 1./sqrt(ae);
N = 2*fix(N/2);  % make N even
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
