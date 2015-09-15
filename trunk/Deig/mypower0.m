function [v,l]=mypower0(A,vn)
% This function returns the largest eigenvalue and eigenvector of A. 
% However, if x0 is given and is the corresponding eigenvector of the
% largest eigenvalue, it will converge to the second largest eigenvalue and
% the corresponding eigenvector. 
%

if nargin >1
   v0 = vn;
   normv02 = norm(v0)^2;
end
n  = size(A,1);
if nargin==1
    vn = rand(n,1);
end

v   = zeros(n,1);
vkn = vn/norm(vn);

it    = 0;
maxit = 1e4;
while(and(norm(vkn-v)>1e-10,it<maxit))
    %norm(vkn-v)>1e-10
   it        = it+1;
   v         = vkn; 
   vk        = A*v;
   if nargin > 1
       vk = vk - (vk'*v0)/normv02*v0;   
   end   
   vkn       = vk/norm(vk);   
end

disp(it)
v = vkn;
l = norm(A*v);