function [var2,M]=diffsimu(n,X,k,nit)

var2    = [];
t       = 10;
dx      = 1/n;
 
xm      = min(min(X));
xM      = max(max(X));

[Xs,Ys] = meshgrid(dx/2:dx:1-dx/2,dx/2:dx:1-dx/2);
Xs      = Xs-0.5;
Ys      = Ys-0.5;



M     = exp(-(((Xs./dx).^2)).*(k*((Ys./dx).^2)));

M     = M/max(max(M));
M     = fftshift(M);

%close all
%figure
%hold on
%surf(M)
%surf(fftshift(M))


for i = 1:nit
   
  %imagesc(X,[xm xM]);
  %axis equal         
  %set(gca,'ydir','normal');
  %plot(X(:,50))
  %axis([1 100 -1 1])
  %Mm(i) = getframe;

  Xf  = fft2(X);
  X   = real(ifft2(Xf.*M));
  var2 = [var2 var(reshape(X,n^2,1))];

end
