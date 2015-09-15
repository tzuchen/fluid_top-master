function [xp,yp] = tentmap2d(x,y,period,i,param)
% functoion for standard map
% period = 1, periodical with period 1
%          0, not periodical
% i      > 0, standard map
%        < 0, inverse standard map
%
%  Tzu-Chen Liang 7-5-2007
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i>0 
   if period ==1
      xp = mod(tentmap(x)+y,1);
      yp = mod(x,1);
   else
      xp = tentmap(x)+y;
      yp = x;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   if period ==1
      xp = mod(y,1);
      yp = mod(x-tentmap(y),1);
   else
      xp = y;
      yp = x-tentmap(y);
   end
end
