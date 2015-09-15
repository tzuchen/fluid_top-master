% demomixingchannelgrid
close all
figure
hold on
box on


fz = 14

load yefilew_save
load zefilew_save

nx = 50
dx = 1/nx;


ym = reshape(yefilew,nx+1,nx+1);
zm = reshape(zefilew,nx+1,nx+1);
[yi,zi] = meshgrid(0:dx:1,0:dx:1);


for i = 1:nx
  for j= 1:nx
   plot([ym(i,j) ym(i+1,j)]  ,[zm(i,j),zm(i+1,j)],'r','linewidth',1.5);
   plot([ym(i,j) ym(i  ,j+1)],[zm(i,j),zm(i  ,j+1)],'r','linewidth',1.5);
  end
end

for i = 1:nx
  for j= 1:nx
   plot([yi(i,j) yi(i+1,j)]  ,[zi(i,j),zi(i+1,j)]);
   plot([yi(i,j) yi(i  ,j+1)],[zi(i,j),zi(i  ,j+1)]);
  end
end

axis equal
axis tight
xlabel('x','fontsize',fz);
ylabel('y','fontsize',fz);
set(gca,'fontsize',fz);
