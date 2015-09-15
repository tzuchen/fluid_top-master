function [xp,yp] = arnoldscatmap(x,y,period,i,param)
% functoion for arnoldscatmap map
% period = 1, periodical with period 1
%          0, not periodical
% i      > 0, arnoldscatmap map
%        < 0, arnoldscatmap standard map
%
%  Tzu-Chen Liang 7-10-2006
%

if nargin > 4
   K = param{1};
else
   K = 0.1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if i>0 
   if period ==1
      xp = mod(2*x + y + K/2/pi*sin(2*pi*x),1);
      yp = mod(  x + y + K/2/pi*sin(2*pi*x),1);
   else
      xp = 2*x + y + K/2/pi*sin(2*pi*x);
      yp =   x + y + K/2/pi*sin(2*pi*x);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   if period ==1
      xp = mod( x-  y,1);
      yp = mod(-x+2*y-K/2/pi*sin(x-y) ,1);
   else
      xp =  x  -y;
      yp = -x+2*y-K/2/pi*sin(x-y);
   end
end
