%demobntncutof
load nlist

close all
figure 
hold on
dlist = [ 100 200 400 1000 2000];
fz = 14;
fzl=18
obj = 1;  %obj=0 paper %obj=1 slides
dlen = length(dlist);

if obj == 0
   cm = zeros(dlen,3);
   ls ={'-','--','-.','d','s'}; 
else
   cm = lines(dlen);
   ls ={'-','-','-','-','-'}; 
end    

for i = 1:dlen
   h(i) = plot(nlist(i,:));
   set(h(i),'color',cm(i,:));
   set(h(i),'linewidth',2);
   set(h(i),'linestyle',ls{i});
end

axis([0 8000 0 1.3])
text(3000,-0.105, '$k$ (iteration)','fontsize',fzl,'interpreter','latex');
ylabel('$|\omega^k_n$-$\bar{\omega}|_{TV}$','fontsize',fzl,'interpreter','latex');
set(gca,'fontsize',fz);
box on


cl = 2000;
cr = 6200;
g(1)= plot([cl cl],[0 1.15],'color',[0 0 0]);
g(2)= plot([cr cr],[0 1.15],'color',[0 0 0]);
%g(3)= plot([cl cr],[1.13 1.13],'color',[0 0 0]);

ct = 3600;
g(3)= plot([ct ct],[0 1.25],'color',[0 0 0]);
%g(5)= plot([0  ct],[1.22 1.22],'color',[0 0 0]);

text(1000,1.25,'cutoff time','interpreter','latex','fontsize',fzl);
text(3800,1.14 ,'decay time','interpreter','latex','fontsize',fzl);
title('Random walk on an $n$-dimensional hypercube','fontsize',fzl,'interpreter','latex')

annotation1 = annotation('doublearrow',[0.329 0.7339],[0.7905 0.7905]);
annotation1 = annotation('doublearrow',[0.138 0.4821],[0.8595 0.8595]);
set(gca,'fontsize',fz)

