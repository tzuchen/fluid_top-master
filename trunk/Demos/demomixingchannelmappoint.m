% demo mixiningchannelmappoint
% 
% Tzu-Chen Liang 10-21-2007

close all
figure
hold on
box on
fz = 14;


load yefilew_save
load zefilew_save

h = plot(yefilew_save,zefilew_save,'.')
set(h,'markersize',15,'linewidth',2);
axis equal
axis tight


xlabel('x','fontsize',fz)
ylabel('y','fontsize',fz)
set(gca,'fontsize',fz);


