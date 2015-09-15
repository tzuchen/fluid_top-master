function [v,w] = lininterp3d(x,vc)
% Given the function values at the 8 corners of a unit cube,
% this function does linear interpolation to find v(x), and 
% it returns the weighting of each corner.  
%
%
% The order of vc is as the following
%  1. (0 0 0)
%  2. (0 0 1)
%  3. (0 1 0)  
%  4. (0 1 1)
%  5. (1 0 0) 
%  6. (1 0 1)
%  7. (1 1 0)
%  8. (1 1 1)
%
% Tzu-Chen Liang  11/15/2005
  
   xn = 1-x;

   c1 = [xn(1,:) ; xn(1,:) ; xn(1,:) ; xn(1,:) ;  x(1,:) ;  x(1,:) ;  x(1,:) ; x(1,:)];
   c2 = [xn(2,:) ; xn(2,:) ;  x(2,:) ;  x(2,:) ; xn(2,:) ; xn(2,:) ;  x(2,:) ; x(2,:)];
   c3 = [xn(3,:) ;  x(3,:) ; xn(3,:) ;  x(3,:) ; xn(3,:) ;  x(3,:) ; xn(3,:) ; x(3,:)];

   w = c1.*c2.*c3;
   v = sum(w.*vc);


