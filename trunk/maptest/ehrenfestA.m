function [A,xe] = ehrenfestA(d)

A = sparse(d+1,d+1);
n = d+1;

for i = 1:d+1
  A(i,i) = 1/(d+1);
end

for i = 2:n
  A(i,i-1) = (i-1)/(d+1);
end

for i = 1:n-1

  A(i,i+1) = (d-(i-1))/(d+1);
end


for i = 0:d
 % xe(i+1) = nchoosek(d,i)/2^d;
end
%xe = xe';


%x0 = zeros(n,1);
%x0(1) = 1;
%x = x0;
%nlist = [];


%for i = 1:1000
%   x = A'*x;  
%end
%xe = x;

%x0 = zeros(n,1);
%x0(1) = 1;
%x = x0;
%nlist = [];
%for i = 1:100
%   x = A'*x;
%   nlist = [nlist sum(abs(x-xe))];
%   %nlist = [nlist sqrt(sum(xe.*(x./xe).^2))];
%  
%end
xe=[];
%plot(nlist)
