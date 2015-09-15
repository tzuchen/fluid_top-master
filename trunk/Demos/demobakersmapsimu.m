%demobakersmapsimu

close all
fzl=24
n=500;
[X,Y] = meshgrid((0:n)/n,(0:n)/n);

Z = cos(2*pi*Y);

imagesc(Z)
axis equal
grid off
box on
axis tight
colormap(gray)
set(gca,'xtick',[])
set(gca,'ytick',[])
text(0,-20,'$k=0$','fontsize',fzl,'interpreter','latex')

figure
Z = cos(2*2*pi*Y);
imagesc(Z)
axis equal
grid off
box on
axis tight
colormap(gray)
set(gca,'xtick',[])
set(gca,'ytick',[])
text(0,-20,'$k=1$','fontsize',fzl,'interpreter','latex')

figure
Z = cos(4*2*pi*Y);
imagesc(Z)
axis equal
grid off
box on
axis tight
colormap(gray)
set(gca,'xtick',[])
set(gca,'ytick',[])
text(0,-20,'$k=2$','fontsize',fzl,'interpreter','latex')
