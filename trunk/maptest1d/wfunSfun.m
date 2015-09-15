function w = wfunSfun(y,k,mun,S,fun,funp) 
%  Calculate w_(mun)^k by given mun and the map S 
%  and given the initial distribution function  
% 
%   0<x<mun   w^0(x) = 1/mun
%   mun<x<1   w^0(x) = 0
%
%  Tzu-Chen Liang 8-2-2007



  if k>0

     yp = feval(S,y,'i');
     wp = 1./abs(feval(S,yp(1,:),'d')).*wfunSfun(yp(1,:),k-1,mun,S,fun,funp);
     wm = 1./abs(feval(S,yp(2,:),'d')).*wfunSfun(yp(2,:),k-1,mun,S,fun,funp);
     w  = (wp+wm);

   w(find(isnan(w))) = 0;

  else       
      w = 0*y;
      ind = find(y<mun);
      w(ind) =feval(funp,y(ind))./feval(fun,mun);

  end
