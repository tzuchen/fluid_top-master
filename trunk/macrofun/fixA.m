function A = fixA(A)
% if the column sum of A is not 1
% this function makes it 1 by averaging neighbor columns!

p = find(sum(A')<0.9);
for i = p
   A(i,:) = (A(i-1,:)+A(i+1,:))/2;
   disp('one column of A is fixed!')
end


 
