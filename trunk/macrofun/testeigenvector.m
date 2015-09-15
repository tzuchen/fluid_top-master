function testeigenvector(A)

OPTS.disp = 0;
[V,D]       = eigs(A',6,'LM',OPTS);
[Dr,I]      = sort(diag(abs(D)));

D = D(I,I);
V = real(V(:,I));

sa = sqrt(size(A,1));

for i = 1:6 
    subplot(3,2,i)
    contourf(reshape(V(:,i),sa,sa));
    title(sprintf('lambda = %f',D(i,i)) )
    axis equal
end
