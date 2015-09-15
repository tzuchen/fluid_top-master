%demo map plot
%
%
%
%  Tzu-Chen Liang 10-7-2007


 n  = 20;
 dx = 1/n;
 dy = dx;
 param{1} = 0.3

 [x0,y0] = meshgrid(-1.5:dx:1.5,-1.5:dy:1.5);
 [xe,ye] = standardmap(x0,y0,0,1,param);
 nd = size(x0,1)-1;


close all
figure 
hold on

%plot(x0,y0,'.b')
%plot(xe,ye,'.r')

for i = 1:nd
   for j = 1:nd
     plot([x0(i,j),x0(i,j+1)],[y0(i,j),y0(i,j+1)])
     plot([x0(i,j),x0(i+1,j)],[y0(i,j),y0(i+1,j)])
  end
end

for i = 1:nd
   for j = 1:nd
   
        plot([xe(i,j),xe(i,j+1)],[ye(i,j),ye(i,j+1)],'r','linewidth',1.5);
        plot([xe(i,j),xe(i+1,j)],[ye(i,j),ye(i+1,j)],'r','linewidth',1.5);
     end
end


box on
axis equal
axis([0 1 0 1])
xlabel('$x_1$','fontsize',16,'interpreter','latex')
ylabel('$x_2$','fontsize',16,'interpreter','latex')
title('Standard map, $\epsilon=0.3$','fontsize',16,'interpreter','latex')
set(gca,'fontsize',16)
