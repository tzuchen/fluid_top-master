load star
C = A;
clear A b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/tzuchen/snopt7111/snopt-dev/matlab
addpath /home/tzuchen/petsc/petsc-2.3.0/bin/matlab
addpath /home/tzuchen/proj/fluid_top/trunk/matlab_send

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = size(C,1);
 
for i = 1:m
  A{i}      = sparse(m,m);
  A{i}(i,i) = 1; 
  b(i)      = 1;
 end
  b = b';

option.parallel = 1; 
 
[X,y,S] = flirtu(C,A,b,option);

[V,D]   = eig(X);
 V      = V*sqrtm(D);

q  = rand(m,1)-0.5;
s  = zeros(m,1);
w  = V*q;
s(find(w>0)) = 1;
s(find(w<0)) = -1;

scatter(x(1,:),x(2,:));
hold on

for i = 1:m
   for j = i+1:m
      if C(i,j)>0
         plot([x(1,i) x(1,j)],[x(2,i) x(2,j)]);
      end
   end
end

h = scatter(x(1,find(s>0)),x(2,find(s>0)),'ro','filled');
set(h,'linewidth',3)
h = scatter(x(1,find(s<0)),x(2,find(s<0)),'go','filled');
set(h,'linewidth',3)
