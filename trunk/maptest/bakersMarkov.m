function A = bakersMarkov(n)
%
% This function generate a Markov matrix for baker's map
% We use this because the map2markov.m doesn't work for 
% baker's map
% n has to be even.
% 
% Tzu-Chen Liang 5-9-2006 



A = sparse(n^2,n^2);

ind = 0;
for j = 1:n
  for i = 1:n
    ind = ind + 1;   
     if i<=n/2;
    %   if mod(j,2)==0
        ij = (floor((j+1)/2)-1)*n+(2*i-1);
        A(ind, [ij ij+1])   = [0.5 0.5];
    %   else
    %    ij = (floor((j+1)/2)-1)*n+(2*i-1);
    %    A(ind, [ij ij+1])   = [0 1];
    %   end
     else 
     %  if mod(j,2)==0
        ij = (floor((j+1)/2)+n/2-1)*n+ (2*i-1-n);
        A(ind, [ij ij+1])   = [0.5 0.5];
     %  else
     %    ij = (floor((j+1)/2)+n/2-1)*n+ (2*i-1-n);
     %    A(ind, [ij ij+1])   = [1 0];
     %  end
     end
   end
end
 
