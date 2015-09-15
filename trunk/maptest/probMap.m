function x = probMap(A,M,p,x)
% 
%  x = A*M^p *x
%
% Tzu-Chen Liang 4-7-2006


    for i = 1:p
       x = M*x;
    end
       x = A*x;


  
