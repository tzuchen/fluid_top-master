function  B = Bgenerate(n,m)
% This function generate the observation matrix B
% 
%
%
% Tzu-Chen Liang 2-6-2007

[X,Y] = meshgrid(1:n,1:n);
blk   = n/m;

X = reshape(ceil(X/blk),n^2,1);
Y = reshape(ceil(Y/blk),n^2,1);

ind = sub2ind([m,m],Y,X);

B = sparse(1:n^2,ind,ones(1,n^2),n^2,m^2);

