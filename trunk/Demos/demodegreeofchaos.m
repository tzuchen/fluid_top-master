% demodegreeofchaos
function []=demodegreeofchaos(x)


close all
x= rand(1e3,1)/10;
y= rand(1e3,1)/10

x1=x;x2=x;x3=x;
y1=y;y2=y;y3=y;

figure

p = 1;

for i = 1:3
 
 [b,l]=ind2sub([3,3],p)
 p = p + 1;
 axes('position',[(l-1)*1/3+0.02  ((b-1)*1/3)+0.02  0.30  0.30])
 h=plot(x3,y3,'.');
 axis equal
 box on
 axis tight
 set(gca,'XTick',[],'XTickLabel',[])
 set(gca,'YTick',[],'YTickLabel',[])
 %set(h,'markersize',1e-0)
 axis([0 1 0 1])
 

 [b,l]=ind2sub([3,3],p)
 p = p + 1;
 axes('position',[(l-1)*1/3+0.02  ((b-1)*1/3)+0.02  0.30  0.30])
 h=plot(x2,y2,'.');
 axis equal
 box on
 axis tight
 set(gca,'XTick',[],'XTickLabel',[])
 set(gca,'YTick',[],'YTickLabel',[])
 %set(h,'markersize',1e-0)
 axis([0 1 0 1])

 [b,l]=ind2sub([3,3],p)
 p = p + 1;
 axes('position',[(l-1)*1/3+0.02  ((b-1)*1/3)+0.02  0.30  0.30])
 h=plot(x1,y1,'.');
 axis equal
 box on
 axis tight
 set(gca,'XTick',[],'XTickLabel',[])
 set(gca,'YTick',[],'YTickLabel',[])
 %set(h,'markersize',1e-0)
 axis([0 1 0 1])

  for j = 1:2
  [x1,y1] = S1(x1,y1);
  [x2,y2] = S2(x2,y2);
  [x3,y3] = S3(x3,y3);
  end
  
end

%title('Ergodic, Mixing and Exact map')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ergodic map
function [xp,yp] = S1(x,y)

xp = mod(sqrt(2)+x,1);
yp = mod(sqrt(3)+y,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mixing map
function [xp,yp] = S2(x,y)

xp = mod(x+y,1);
yp = mod(x+2*y,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact
function [xp,yp] = S3(x,y)

xp = mod(3*x+y,1);
yp = mod(x+3*y,1);


