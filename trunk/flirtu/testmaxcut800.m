
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath /home/tzuchen/snopt7111/snopt-dev/matlab
addpath /home/tzuchen/petsc/petsc-2.3.0/bin/matlab
addpath /home/tzuchen/proj/fluid_top/trunk/matlab_send

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxcut800
m = size(C,1);
 


option.parallel = 0; 

m = 800
C = C(1:m,1:m);
clear A
b = b(1:m);
for i=1:m
C(i,i) = 1;
A{i} = sparse(m,m);
A{i}(i,i)=1;
end

 
[X,y,S] = flirtu(C,A,b,option);

