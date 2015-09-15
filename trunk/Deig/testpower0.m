% testpower0
% ASSUMPTION: the largest 2 eigenvalues are real and positive!
n = 1000;
A = rand(n,n);
B = diag(sparse(ones(n-1,1)),1)+diag(sparse(ones(n-2,1)),2)+diag(sparse(ones(n-3,1)),3);
B = sparse(B+B' + speye(n));
A = sparse(A.*B);

for i= 1:n
 A(i,:)=A(i,:)/sum(A(i,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the largest 2 eigenvalues and eigenvectors
[V,D]   = eigs(A);
d       = diag(D);
[nouse,Ind] = sort(d);

maxeig1 = d(Ind(end))
maxeig2 = d(Ind(end-1))
maxv1   = V(:,Ind(end));
maxv2   = V(:,Ind(end-1));

disp(sprintf('The maximum 2 eigenvalues are %d and %d',maxeig1,maxeig2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[v1,l1]=mypower0(A);

a1 = A(1,:);
A2 = A- sparse(v1*a1*l1/(a1*v1));

[v2,l2]=mypower0(A2);


disp(maxeig1)
disp(l1)
disp(maxeig2)
disp(l2)