%demo normcompare plot
% Tzu-Chen Liang 10-18-2007

   fz = 14;
   fzl = 16;
   obj = 0;  %obj=1 slide, obj=0 paper
   if obj == 1
     cm = jet(6);
     linestylelist= {'-','-','-','-','-','-'};
   else
     cm = zeros(6,3);
     linestylelist= {'-','--','-.','d','x','s'};
   end

load alln
%[Astruc] = maprefine2(500,[],@standardmap,0.3);
%[Xa1,var2,var1,mnorm,snorm,sv,Mm,alln] = testAf(Astruc,[],'cosx',31,0,1,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure
hold on
grid on
box on

p = [1 3 2 4 5];
for i = p

   h(i) = plot(0:30,alln(i,:))
   set(h(i),'linewidth',2,'linestyle',linestylelist{i},'color',cm(i,:),'markersize',6)

end

%legend('$L_1$','$L_2$','mix-norm','H$_{0.5}$','$H_{-1}$','fontsize',fzl)
legend('L_1','L_2','mix-norm','H_{-0.5}','H_{-1}','fontsize',fz)
xlabel('iteration','interpreter','latex','fontsize',fzl)
ylabel('$||.||$, normalized','interpreter','latex','fontsize',fzl)
set(gca,'fontsize',fz)
for i = 4:5
   g(i) = plot(0:30,alln(i,:),'linewidth',1,'color',cm(i,:))
end
title('Standard Map with $\epsilon=0.3$ simulated by $A_n$, $n=500^2$','interpreter','latex','fontsize',fzl)
