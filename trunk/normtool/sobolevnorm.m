function [w] = sobolevnorm(X,domain,p,slist)
%
% Calculate Sobolev norm W(s,p) of X
% where s can be positive, zero or negative
%  
% 
%
%
%   Tzu-Chen Liang 7-7-2006

if prod(size(X)) ==1
  % only calculate the cofficient D;
  n = X;
  P = 4*pi^2*(0:n).^2; 
  w = zeros(length(slist),n+1);
else

   if domain =='t'
     Xf = fft2(X);
   elseif domain == 'f'
     Xf = X;
   else
     error('domain must be t or f');
   end
  n     = size(Xf,1);
  nlist =  [n/2:-1:1,0:1:n/2-1];
  P     = 4*pi^2*(ones(n,1)*nlist.^2+ (ones(n,1)*nlist.^2)');
end



for i = 1:length(slist) 
   s = slist(i);

   if p ==2
      if     s==0   % W(0,2) = H(0) = L(2) 
          D{i} = 1; 
      elseif s > 0  % W(s,2) = H(s)
          D{i} = (1+ P).^s;   
      elseif s < 0  % W(s,2) = H(s)
          D{i} = 1./(1 + P).^(-s);
      else
          error('k must be a real number');
      end
   
   else % p~=2 
       error('p~=2 has not been implemented') 
   end         
end


for i = 1:length(slist) 
    if prod(size(X)) ==1     
         if size(D{i}) ==1
            w(i,:) = D{i}.*ones(1,n+1);
         else
            w(i,:) = D{i};

         end
    else
         D{i} = fftshift(D{i});
         w(i) = sqrt(sum(sum(D{i}.*(abs(Xf).^2))))/n/n;
    end
end











