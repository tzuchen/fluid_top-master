function testsigma(A,x0,V,D)
%close all
nit = 400;
nit = 400;
close all
sa = sqrt(size(A,1));
Xa = zeros(sa,sa);

% Only give A
if nargin == 1
    Xa(1:fix(sa/2),:) = 0;
    Xa(fix(sa/2):end,:) = 1;
    xa0 = reshape(Xa,(sa)*(sa),1);  
else

% Initial distribution is given
    if prod(size(x0))>1;     
       xa0 = x0;
    end

% Use eigen vectors
    if prod(size(x0))==1

      if nargin < 4
         OPTS.disp=0;
         [V,D]       = eigs(A,6,'LM',OPTS);
         [Dr,I]      = sort(diag(abs(D)));
         D = D(I,I);
         V = real(V(:,I));
      end
 
        % use the x0-st eigenvector as initial distribution
     
             v = real(V(:,x0));
      
        % scale it to the positive orthant.
        xa0 = (v-min(v))/(max(v)-min(v));
    
    end


end
  OPTS.disp=0;
         [V,D]       = eigs(A,6,'LM',OPTS);
         [Dr,I]      = sort(diag(abs(D)));

  [dr,I]      = sort(abs(diag(D)));

  xa1 = xa0;




dy = 1/sa;
y  = dy/2:dy:1-dy/2;
z  = dy/2:dy:1-dy/2;
[Y,Z] = meshgrid(y,z);

ds = 4/sa
yi  = ds/2:ds:1-ds/2;
zi  = ds/2:ds:1-ds/2;
[Yi,Zi] = meshgrid(yi,zi);

sig=[];
       
   for i = 1:nit
 
        xa1     = A*xa1;          
        Xa1    = reshape(xa1,sa,sa);
        Xas    = interp2(Y,Z,Xa1,Yi,Zi);
        xas    = reshape(Xas,prod(size(Xas)),1);
        sig(i) = std(xas);

    end
  

drlist = zeros(length(dr),nit);

for i = 1:nit
   dr2(i,:) = sig(1)*dr.^(i-1);
end



h = semilogy(sig)
set(h,'linewidth',2);
hold on
   semilogy(dr2)
grid on
xlabel('Number of Iterations')
ylabel('\sigma')






     
