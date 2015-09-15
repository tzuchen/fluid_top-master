function [xp,yp] = logisticmap2d(x,y,period,i,param)
% functoion for standard map
% period = 1, periodical with period 1
%          0, not periodical
% i      > 0, standard map
%        < 0, inverse standard map
%
%  Tzu-Chen Liang 7-5-2007
%

if nargin > 4
   c1 = param{1};
else
   c1 = 0.1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i>0 
   if period ==1
      xp = mod(4*x.*(1-x)+y,1);
      yp = mod(x,1);
   else
      xp = 4*x.*(1-x)+y;
      yp = x;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   if period ==1
      xp = mod(y,1);
      yp = mod(x-4*y.*(1-y),1);
   else
      xp = y;
      yp = x-4*y.*(1-y);
   end
end
