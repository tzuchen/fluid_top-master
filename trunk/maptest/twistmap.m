function [xf,yf] = twistmap(x,y,period,i,param)

xc = param{1};
yc = param{2};
dir = param{3};

xspec = 0:0.1:1;
yspec = 0:0.1:1;
xmspec = interp1([0 0.5 1],[0 xc 1],xspec,'cubic');
ymspec = interp1([0 0.5 1],[0 yc 1],yspec,'cubic');


  
  
  x = interp1(xmspec,xspec,x,'linear');
  y = interp1(ymspec,yspec,y,'linear');

[th,r] = cart2pol(x-0.5,y-0.5);
 rs    = rcfind(0.5,0.5,th);

dx  = (x-x.^2).*(1-2*y);
dy  = (y-y.^2).*(1-2*x);
dv  = (dx.^2+dy.^2).^0.5+1e-10;
dx  = dx./dv;
dy  = dy./dv;
t   = 0.1*(r./rs.*(1-r./rs));
%t   = 0.1*(r./rs.*(1-r./rs)).^2;
xf  = x+dir*t.*dx;
yf  = y-dir*t.*dy;

dxf  = (xf-xf.^2).*(1-2*yf);
dyf  = (yf-yf.^2).*(1-2*xf);
dvf  = (dxf.^2+dyf.^2).^0.5+1e-10;
dxf  = dxf./dv;
dyf  = dyf./dv;

xf  = x+dir*t.*dx;
yf  = y-dir*t.*dyf;


 

  xf(find(xf<0))=0;
  yf(find(yf<0))=0;
  xf(find(xf>1))=1;
  yf(find(yf>1))=1;
   



xf = interp1(xspec,xmspec,xf,'linear');
yf = interp1(yspec,ymspec,yf,'linear');
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = rcfind(xc,yc,th)

  r1 = abs((0.5+sign(cos(th))*(0.5-xc))./cos(th));
  r2 = abs((0.5+sign(sin(th))*(0.5-yc))./cos(pi/2-th));

  r = min(r1,r2);
  

  
  
    
        
