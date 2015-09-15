function [xp,yp] = tentmap2d2(x,y,period,i,param)
% functoion for tentmap2d2 map
% period = 1, periodical with period 1
%          0, not periodical
% i      > 0, standard map
%        < 0, inverse standard map
%
%  Tzu-Chen Liang 7-5-2007
%

 xp = x;
 yp = y;
 

if i > 0
   ind1 = find(and(0<=x,x<1/2));
   ind2 = setdiff(1:prod(size(x)),ind1);
    xp(ind1)  = 2*x(ind1);
    yp(ind1)  = y(ind1)/2;
    xp(ind2)  = 2-2*x(ind1);
    yp(ind2)  = 1-y(ind2)/2;%(y(ind2)+1)/2;
else
   ind1 =find(and(0<=y,y<1/2));
   ind2 = setdiff(1:prod(size(x)),ind1);
    xp(ind1) = x(ind1)/2;
    yp(ind1) = 2*y(ind1);
    xp(ind2) = 1-x(ind2)/2;
    yp(ind2) = 2-2*y(ind2);%2*y(ind2)-1;     
end

 if period ==1
   xp  = mod(xp,1);
   yp  = mod(yp,1);
 end
