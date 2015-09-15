function [A,b]=pt2ab(x,y);
% given points x,y in 2-D (CW order), this function returns 
% A and b such that Ax < b represents the polygon formed by
% these points. 
% Tzu-Chen Liang 1-19-2006


n    = length(x)+1;
x(n) = x(1);
y(n) = y(1);
A    = zeros(n-1,2);
b    = zeros(n-1,1);

for i = 1: n-1

   if y(i+1)==y(i)
     A(i,1) = 0;
     A(i,2) = sign(x(i+1)-x(i));
     b(i)   = sign(A(i,2))*y(i);
   else

     A(i,1) = -sign(y(i+1)-y(i));
     A(i,2) = -(x(i+1)-x(i))/(y(i+1)-y(i))*sign(A(i,1)); 
     b(i)      = (A(i,1)*x(i)+A(i,2)*y(i));
   end


end

