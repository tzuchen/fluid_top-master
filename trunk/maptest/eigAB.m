% given A1,A2
% calculate the invariant distribution

n = size(A1.A,1);
p = 100

x = rand(n,1);

for i= 1:5000
  xp = x;
  if mod(i,p)>p/2
   x = A1.A'*x;
  else
   x = A2.A'*x;
 end
  norm(x-xp)
end

