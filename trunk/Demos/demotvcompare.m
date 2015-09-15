%test

eps = 0.5
r   = 0.9

p = 3

d(1) = 0.5+eps;
for i = 2:p
  d(i) = d(i-1)*r;  
end

 dn = 1-d; 

 wmat = zeros(2^p,p);
 k = 1;
 
for i = 1:2^p
  for j = 1:p
  
     if sin((i-0.5)/2^p*2*pi*k) >0
        wmat(i,j) = d(j);
     else
        wmat(i,j) = dn(j);
     end
     
  end
k = k*2;
end


