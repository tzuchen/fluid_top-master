function w = wfun(y,k,mun) 
%  Calculate w_(mun)^k by given mun 
%   
%   0<x<mun   w^0(x) = 1/mun
%   mun<x<1   w^0(x) = 0
%
%  Tzu-Chen Liang 5-9-2007





  if k>0

     yp = (1-(1-y).^0.5)./2;
     ym = (1+(1-y).^0.5)./2;

     wp = abs(1./(4-8*yp)).*wfun(yp,k-1,mun);
     wm = abs(1./(4-8*ym)).*wfun(ym,k-1,mun);
     w  = (wp+wm);

  else       
      w = 0*y;
      ind = find(y<mun);
      w(ind) =1/mun;
  end


