fz = 14;
fzl = 16;
obj = 0;  %obj=1 slide, obj=0 paper
close all

if obj == 1
 cm =lines(6);%jet(6);
 linestylelist= {'-','-','-','-','-','-'};
else
 cm = zeros(6,3);
 linestylelist= {'-','--','-.','d','s','x'};
end



close
figure
hold on
box on
grid on
k = 0:50
Dlist = [1e-6 1e-8 1e-10 1e-12 1e-14]
clear nlist

for i=1:length(Dlist)
D = Dlist(i);
nlist(i,:)=exp(-4*pi^2*D*2.^(2*k)).^2;
h = plot(0:length(k)-1,nlist(i,:));
 set(h,'linewidth',2,'linestyle',linestylelist{i},'color',cm(i,:))
end

for i=4:length(Dlist)
  hp = plot(0:length(k)-1,nlist(i,:));
   set(hp,'linewidth',1,'color',[0,0,0])
end


for i = 1:length(Dlist)
   mylegend{i} = sprintf('D = %1.1e', Dlist(i));
end
axis([0 30 0 1.2])
legend(mylegend,'interpreter','latex')

title('Baker''s Map Cutoff','fontsize',16,'interpreter','latex')
xlabel('$k$ (iteration)','fontsize',16,'interpreter','latex')
ylabel('var($f$)','fontsize',16,'interpreter','latex')
set(gca,'fontsize',14)

iter = 51;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
noval = 0.5

for i = 1:dl
  for j = 1:iter-1
    if nlist(i,j+1)>=nlist(i,j);
      nlist(i,j+1) = nlist(i,j)-1e-9; 
    end
  end
end


for i = 1:dl
  itern(i)   = interp1(nlist(i,:),0:iter-1,noval);
end



figure
hold on

for i = 1:dl
  h(i) = plot([0:iter-1]/itern(i),nlist(i,:));
  set(h(i),'color',cm(i,:));  
end




legend(w,'fontsize',fz,'interpreter','latex');
axis([0 2 0 1.2])
set(h,'linewidth',2,'markersize',6);

box on
grid on
xlabel('k_n/t_n normalized iteration','fontsize',fzl,'interpreter','latex')
ylabel('var($f$)','fontsize',fzl,'interpreter','latex')
title('Baker''s Map Cutoff','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',fz);

