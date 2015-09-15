function [xp,yp] = standardmap(x,y,period,i,param)
% functoion for standard map
% period = 1, periodical with period 1
%          0, not periodical
% i      > 0, standard map
%        < 0, inverse standard map
%
%  Tzu-Chen Liang 3-24-2006
%

if nargin > 4
   c1 = param{1};
else
   c1 = 0.1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i>0 
   if period ==1
      xp = mod(x + y + c1*sin(2*pi*x),1);
      yp = mod(y+ c1*sin(2*pi*x),1);
   else
      xp = x + y + c1*sin(2*pi*x);
      yp = y+ c1*sin(2*pi*x);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   if period ==1
      xp = mod(x-y,1);
      yp = mod(y-c1*sin(2*pi*(x-y)),1);
   else
      xp = x-y;
      yp = y-c1*sin(2*pi*(x-y));
   end
end
