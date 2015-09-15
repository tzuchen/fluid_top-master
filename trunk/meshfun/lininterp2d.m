function [v,w] = lininterp2d(x,vc)
% Given the function values at the 4 corners of a unit square,
% this function does linear interpolation to find v(x), and 
% it returns the weighting of each corner.  
%
%
% The order of vc is as the following
%  1. (0 0)
%  2. (0 1)
%  3. (1 0)  
%  4. (1 1)
%
% Tzu-Chen Liang  11/16/2005
  
   xn = 1-x;

   c1 = [xn(1,:) ; xn(1,:) ;  x(1,:) ;  x(1,:)];
   c2 = [xn(2,:) ;  x(2,:) ; xn(2,:) ;  x(2,:)];
   
   w = c1.*c2;
   v = sum(w.*vc);
