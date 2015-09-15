function [xf,yf] = twistmap2(x,y,period,i,param)

xc  = param{1};
yc  = param{2};
r0  = param{3};
dir = param{4};


[th,r]= cart2pol(x-xc,y-yc);
ind = find(r<r0);

v   = 500*((r-r0).*r).^2;


th(ind)  = th(ind) + dir*v(ind);

[xs,ys] = pol2cart(th,r);


xf = xs + xc;
yf = ys + yc;

 

  xf(find(xf<0))=0;
  yf(find(yf<0))=0;
  xf(find(xf>1))=1;
  yf(find(yf>1))=1;
   



