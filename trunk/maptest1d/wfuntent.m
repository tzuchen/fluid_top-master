function w = wfuntent(y,k,mun) 
%  Calculate w_(mun)^k by given mun 
%   
%   0<x<mun   w^0(x) = 1/mun
%   mun<x<1   w^0(x) = 0
%
%  Tzu-Chen Liang 5-9-2007


  if k>0

     yp = 1-y/2;
     ym = y/2;

     wp = 1./2.*wfuntent(yp,k-1,mun);
     wm = 1./2.*wfuntent(ym,k-1,mun);
     w  = (wp+wm);

  else       
      w = 0*y;
      ind = find(y<mun);
       w(ind) =1/mun;
      %w = pi.*sin(pi*y)/2;
      %w = w/sin(mun*pi/2)^2;
      %w(find(y>mun)) = 0;
  end
