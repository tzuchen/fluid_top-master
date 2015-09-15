function  [Dlist,ld]= Dvslambda(x0,P,n,statelist,Dstar);
%
% 
%
%
% Tzu-Chen Liang  8-1-2006

Dmax = 1e-2%100*Dstar;
Dint = (Dmax - Dstar)/10;
Dlist = Dstar:Dint:Dmax;


figure 
hold on
for i = 1:length(Dlist)
  Pd = freqmapaddD(P,n,statelist,Dstar,Dlist(i));
  x = x0;  
  lognormx = [];
   for k = 1:200
     x = Pd*x;
     lognormx = [lognormx log(norm(x))];
   end
    
     ld(i)= -(lognormx(end)-lognormx(end-10+1))/10;

end





