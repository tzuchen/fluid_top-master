function xl = testmap(S)

x = rand

xl = zeros(10000,1);
for i = 1:10000
  xl(i) = x;
  x = feval(S,x);

end

