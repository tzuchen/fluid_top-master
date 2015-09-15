function [M]= mapsimu(n)

dx = 1/n;
dy = 1/n;

x = dx/2:dx:1-dx/2;
y = dy/2:dy:1-dy/2;

[X,Y]           = meshgrid(x,y);
xloc   = reshape(X,n^2,1);
yloc   = reshape(Y,n^2,1);

rind  = find(xloc>0.5);
bind  = setdiff(1:n^2,rind);

xrloc = xloc(rind);
yrloc = yloc(rind);
xbloc = xloc(bind);
ybloc = yloc(bind);
 


close all
figure
hold on



for i = 1:100 
   
   %plot(xrloc,yrloc,'.r')
   scatter(xrloc,yrloc,50,'r','filled')
   hold on
   %plot(xbloc,ybloc,'.b')
   scatter(xbloc,ybloc,50,'b','filled')
   [xrloc,yrloc] = standardmap(xrloc,yrloc); 
   [xbloc,ybloc] = standardmap(xbloc,ybloc);
   hold off
   M(i) = getframe;
   
end
movie(M)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y] = henonmap(x,y)

c1 = 0.1;
c2 = 1; 
x = mod(y,1);
y = mod(-c1*x +c2 - y.^2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,y] = standardmap(x,y)

c1 = 0.04;

x = mod(x-y,1);
y = mod(y-c1*sin(2*pi*(x-y)),1);



