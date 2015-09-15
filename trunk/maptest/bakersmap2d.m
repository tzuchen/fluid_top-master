function [xp,yp] = bakersmap2d(x,y,period,i,param)
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
      xp = mod(S(x)+y,1);
      yp = mod(x,1);
   else
      xp = S(x)+y;
      yp = x;
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   if period ==1
      xp = mod(y,1);
      yp = mod(x-S(y),1);
   else
      xp = y;
      yp = x-S(y);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=S(x)
y = mod(2*x,1);
