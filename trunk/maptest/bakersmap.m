function [xp,yp] = bakersmap(x,y,period,i,param)
% functoion for baker's map
% period = 1, periodical with period 1
%          0, not periodical
% i      > 0, baker's map
%        < 0, inverse baker's map
%
%  Tzu-Chen Liang 5-9-2006
%
 
 xp = x;
 yp = y;
 

if i > 0
   ind1 = find(and(0<=x,x<1/2));
   ind2 = setdiff(1:prod(size(x)),ind1);
    xp(ind1)  = 2*x(ind1);
    yp(ind1)  = y(ind1)/2;
    xp(ind2)  = 2*x(ind1)-1;
    yp(ind2)  = (y(ind2)+1)/2;
else
   ind1 =find(and(0<=y,y<1/2));
   ind2 = setdiff(1:prod(size(x)),ind1);
    xp(ind1) = x(ind1)/2;
    yp(ind1) = 2*y(ind1);
    xp(ind2) = x(ind2)/2 +1/2;
    yp(ind2) = 2*y(ind2)-1;     
end

 if period ==1
   xp  = mod(xp,1);
   yp  = mod(yp,1);
 end

  
