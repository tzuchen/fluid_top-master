function [X,A]=RiffleShuffleA(n)

A = zeros(n,n);
a = 2;

A(1,1) =  EulerianNumber(n,1)*(n+1)/2^n;
A(2,1) =  EulerianNumber(n,2)/2^n;



x0 = zeros(n,1);
x0(1) =  1;
X = zeros(n,n+1);
X(:,1) = x0;
for i = 1:n+1 % number of iter
  for j = 1:n % number of rising seq
      a = 2^(i);
      if n+a-j>=n
   X(j,i) =EulerianNumber(n,j)* nchoosek(n+a-j,n)/a^n;

    end
  end
end



for i=1:n-1

    A(:,i+1)= (X(:,i+1)-A*X(:,i))./X(i+1,i)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function knj = EulerianNumber(n,k)
k=k-1;
knj=0;
for j = 0:k+1
   knj = knj + (-1)^j*nchoosek(n+1,j)*(k-j+1)^n;
end
