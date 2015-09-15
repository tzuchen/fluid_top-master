function [normx,A,x] = Asimu(A,iter)
close all
n = size(A,1);

[V,D] = eig(A');
d = abs(diag(D));
[maxd,maxind] = max(d);
v = abs(V(:,maxind));
v = v/sum(v);
x = zeros(n,1);
x(1) = 1;

normx = sqrt(sum(abs(x-v).^2))/n;
normx = (sum(abs(x-v)))/n;

for i =1:iter

  x = A'*x;
  %normx = [normx sqrt(sum(abs(x-v).^2))/n];
  normx = [normx sum(abs(x-v))/n];
end

semilogy(0:iter,normx)
figure
plot(0:iter,normx/normx(1))

% Asimu

x
