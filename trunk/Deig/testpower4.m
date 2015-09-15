% test for my power method

% The size of the matrix
n  = 5;
% the 1-d entries that I want their derivatives.
plist = [1:5];

% a small increment for finite difference method
d  = 0.001; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here we generate a matrix A
A  = rand(n,n);
%A  = A + A';
% 
%  A = [0.4  0.0 0.6;
%       0.2  0.6 0.2;
%       0.4  0.4 0.2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the largest eigenvalue and the eivenvector
[V,D]   = eig(A);
[lanmax,Ind] = max(max(D));
v      = V(:,Ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dv = [];
Dlan = [];
for k = 1:length(plist)
    Ad = A;    
    [i,j] = find(reshape(1:n^2,n,n)==plist(k));
    %[i,j]
    Ad(i,j) = A(i,j) + d;
    [Vd,Dd] = eig(Ad);
    lanmaxd = Dd(Ind,Ind);
    vd       = Vd(:,Ind);

    Dlan = [Dlan (lanmaxd-lanmax)/d];
    Dv   = [Dv (vd-v)/d];
end

    
    [v1,lan1,Dvpower,Dlanpower] = mypower4(A,plist);

dl = lanmax-lan1;
disp(sprintf('eigenvalue difference = %d',dl))
disp(sprintf('The above number indicates the power method converges to the right eigenvalue or not. '))

disp(sprintf('Dv/Daij by difference'))
disp(Dv)
disp(sprintf('Dv/Daij by modified power method'))
disp(Dvpower)

disp(sprintf('Dlambda/Daij by difference'))
disp(Dlan)
disp(sprintf('Dlambda/Daij by modified power method'))
disp(Dlanpower)