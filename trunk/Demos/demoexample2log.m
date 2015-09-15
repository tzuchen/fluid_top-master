close 
figure 
hold on
box on
grid on
for i=2:4

  h(i) =plot(tlist,mydata{i})
  set(h(i),'color',cm(i,:) );
  set(h(i),'linestyle',linestylelist{i},'linewidth',2);

  if obj==1
    ht(i) =plot(tlist,mydata{i},'o')
    set(ht(i),'color',cm(i,:),'markersize',6);
    set(ht(i),'linestyle',linestylelist{i},'linewidth',2);
  end

%h3d(i) =plot(t3dlist,m3d{i})

end





axis([0 3 0 0.55 ])
xlabel('$x$(cm)','fontsize',fzl,'interpreter','latex')
ylabel('Standard Deviation','fontsize',fzl,'interpreter','latex')
title('Mixing trajectories in log scale','fontsize',fzl,'interpreter','latex')

for i = 1:6
  casen{i} = sprintf('Pe = %2.2e',Pelist(i));
end

%legend(h,casen)

for i = 4:6
 % gh(i) =plot(tlist,mydata{i},'linewidth',1,'color',cm(i,:));
end


plot([0 3],[ 0.05 0.05],'linewidth',2,'linestyle','--')
text(0.8,0.066,'$x_{90}$','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',fz)
set(gca,'yscale','log')

