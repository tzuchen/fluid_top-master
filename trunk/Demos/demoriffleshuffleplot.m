%demoriffleshuffleplot
close all
figure
box on

fzl = 16

x = [1 1 1 1 0.924 0.614 0.334 0.167 0.085 0.043];
h = plot(1:10,x);
set(h,'linewidth',2)
hold on
h = plot(1:10,x,'s');
set(h,'linewidth',2,'markersize',10)
xlabel('number of shuffles','fontsize',fzl,'interpreter','latex')
ylabel('How random?','fontsize',fzl,'interpreter','latex')
title('The total variation distance for Riffle Shuffle of $52$ cards','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',14 )
grid on
axis([1 10 0 1.03])
