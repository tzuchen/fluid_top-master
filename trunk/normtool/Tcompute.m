function [T,k] = Tcompute(n);


x = -n:n;
y = -n:n;

[k1,k2] = meshgrid(x,y);

k = [k1(1:end);k2(1:end)];


[ksort,ind]=sort(k(1,:).^2+k(2,:).^2);

k(1,:) = k(1,ind);
k(2,:) = k(2,ind);

s=length(ind);

T = zeros(s,s);

for i= 1:s
  for j = 1:s
   T(i,j) = Tkq(k(:,i),k(:,j));
  end
end
