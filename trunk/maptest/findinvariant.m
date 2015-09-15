function [x] = findinvariant(n,dir,mapfunction,varargin)
%
% This function return an approximate invariant distribution of 
% the given map. the size of x is n*n
%
%
% Tzu-Chen Liang   06-06-2006

 
 k  = 2;
 nk = n*k;
 x  = 1/nk/2:1/nk:(1-1/nk/2);
 y  = 1/nk/2:1/nk:(1-1/nk/2);
 [X,Y] = meshgrid(x,y);
 Xp = mod(reshape(fix(X*n),nk^2,1),n)+1;
 Yp = mod(reshape(fix(Y*n),nk^2,1),n)+1;
 Xvec = sub2ind([n,n],Xp,Yp);

 x = zeros(1,n^2);
 y = x;
 maxi = 400;
 err  = inf;
 param    = varargin;

 figure 
 hold on

 for i = 1:maxi 

  [X,Y]  = feval(mapfunction,X,Y,1,dir,param);


  h = plot(X,Y,'.b');
  set(h,'markersize',1)
  Xp = reshape(floor(X*n),nk^2,1)+1;
  Yp = reshape(floor(Y*n),nk^2,1)+1;
  Xpvec = sub2ind([n,n],Xp,Yp);
  
  
  Xvec  = Xpvec;
 
  if i>10
     x    = x+  hist(Xvec,1:n^2);
     yp   = x/sum(x); 
     err  = norm(yp-y);
     y    = yp;
  end
  
  disp(err)
  if err<1e-5
     break;
  end
  
 end
  x = x/sum(x);

  box on
  axis equal
  axis([0 1 0 1])
  %imagesc(reshape(x,n,n))
  %set(gca,'ydir','normal');

  %figure 
  %h = plot(X,Y,'.b');
  %set(h,'markersize',1)

  

  



