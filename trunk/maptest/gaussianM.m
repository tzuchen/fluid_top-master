function [M,D] = gaussianM(n,sigma)
%
%  This function generates a (n^2)*(n^2) matrix
%  
%
%  Tzu-Chen Liang 4-6-2006

xl = 1;

dx  = xl/n;

x = dx/2:dx:xl-dx/2;
y = dx/2:dx:xl-dx/2;

[X,Y] = meshgrid(x,y);
Xv    = reshape(X,n^2,1);
Yv    = reshape(Y,n^2,1);


M = sparse(n^2,n^2);
D = sparse(n^2,n^2);


ind = 1;
 for j = 1:n
   for i = 1:n
        Dxi = abs(Xv-X(i,j));
        Dyi = abs(Yv-Y(i,j));
        
        D(ind,:) = sqrt((min(Dxi,xl-Dxi)).^2+(min(Dyi,xl-Dyi)).^2);
        mi    = 1/sigma/sqrt(2*pi)*exp(-(D(ind,:).^2)./(2*(sigma^2)));
        nzi   = find(mi>0.1);
        M(ind,nzi) = mi(nzi)/sum(mi(nzi));
        ind = ind+1;
   end
 end

