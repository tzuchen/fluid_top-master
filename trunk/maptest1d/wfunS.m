function w = wfunS(y,k,mun,S) 
%  Calculate w_(mun)^k by given mun and the map S 
%   
%   0<x<mun   w^0(x) = 1/mun
%   mun<x<1   w^0(x) = 0
%
%  Tzu-Chen Liang 5-21-2007



  if k>0

     yp = feval(S,y,'i');
     wp = 1./abs(feval(S,yp(1,:),'d')).*wfunS(yp(1,:),k-1,mun,S);
     wm = 1./abs(feval(S,yp(2,:),'d')).*wfunS(yp(2,:),k-1,mun,S);
     w  = (wp+wm);

   w(find(isnan(w))) = 0;

  else       
      w = 0*y;
      ind = find(y<mun);
      w(ind) =1/mun;
  end

