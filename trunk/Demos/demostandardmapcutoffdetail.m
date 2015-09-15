% demostandardmapcutoffdetail
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fz = 14;
fzl =16;
obj = 0;
line2 = 1;


if obj ==0
cl = [0 0 0;0 0 0];
else
cl = [0 0 1; 1 0 0];
end

close all
figure 
hold on
box on
grid on

Dlist1 = 1./dlist1.^2;
Dlist2 = 1./dlist2.^2;

if line2==1

cf(1) = plot(Dlist1,itern1)
set(cf(1),'linewidth',2,'color',cl(1,:));
df(1) = plot(Dlist1,itern1,'s')
set(df(1),'markersize',10,'linewidth',2,'color',cl(1,:));

cf(2) = plot(Dlist2,itern2,'color',cl(2,:))
set(cf(2),'linewidth',2);
df(2) = plot(Dlist2,itern2,'d','color',cl(2,:))
set(df(2),'markersize',10,'linewidth',2);

set(gca,'xscale','log')
set(gca,'fontsize',fz)
%xlabel('$1/n^2$','fontsize',fzl,'interpreter','latex')
text(10^-8.5,5.25,'$1/n$','fontsize',fzl,'interpreter','latex')
ylabel('$t_n$ (cutoff time)','fontsize',fzl,'interpreter','latex')
axis([10^(-10) 10^(-6.5) 6 14])
legend(df,'Numerical Diffusion','Smoothing Step')
title('Cutoff time versus $1/n$','fontsize',fzl,'interpreter','latex')

else

cf(1) = plot(1./Dlist1,itern1)
set(cf(1),'linewidth',2,'color',cl(1,:));
df(1) = plot(1./Dlist1,itern1,'s')
set(df(1),'markersize',10,'linewidth',2,'color',cl(1,:));
set(gca,'xscale','log')
set(gca,'fontsize',fz)
xlabel('$n$','fontsize',fzl,'interpreter','latex')
text(10^8.5,5.25,'$n$','fontsize',fzl,'interpreter','latex')
ylabel('$t_n$ (cutoff time)','fontsize',fzl,'interpreter','latex')
axis([10^(6.5) 10^(10)  8 14])
%legend(df,'Numerical Diffusion','Smoothing Step')
title('Cutoff time versus $n$','fontsize',fzl,'interpreter','latex')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
hold on
box on
grid on

if line2 ==1

gf(1) = plot(Dlist1,totalarea1-areai1,'color',cl(1,:))
set(gf(1),'linewidth',2);
bf(1) = plot(Dlist1,totalarea1-areai1,'s','color',cl(1,:))
set(bf(1),'markersize',10,'linewidth',2);

gf(2) = plot(Dlist2,totalarea2-areai2,'color',cl(2,:))
set(gf(2),'linewidth',2);
bf(2) = plot(Dlist2,totalarea2-areai2,'d','color',cl(2,:))
set(bf(2),'markersize',10,'linewidth',2);
set(gca,'fontsize',fz);
set(gca,'xscale','log')


%xlabel('$1/n^2$','fontsize',fzl,'interpreter','latex')
text(10^-8.5,0.144,'$1/n$','fontsize',fzl,'interpreter','latex')
ylabel('$\Delta$','fontsize',fzl,'interpreter','latex')
axis([10^(-10) 10^(-6.5) 0.15 0.21])
legend(bf,'Numerical Diffusion','Smoothing Step','location','northwest')
title('$\Delta$ versus $1/n$','fontsize',fzl,'interpreter','latex')

else
  gf(1) = plot(Dlist1,totalarea1-areai1,'color',cl(1,:))
  set(gf(1),'linewidth',2);
  bf(1) = plot(Dlist1,totalarea1-areai1,'s','color',cl(1,:))
  set(bf(1),'markersize',10,'linewidth',2);
  set(gca,'fontsize',fz);
  set(gca,'xscale','log')
  text(10^-8.5,0.144,'$1/n$','fontsize',fzl,'interpreter','latex')
  ylabel('$\Delta$','fontsize',fzl,'interpreter','latex')
  axis([10^(-10) 10^(-6.5) 0.15 0.21])
  %legend(bf,'Numerical Diffusion','Smoothing Step','location','northwest')
  title('$\Delta$ versus $1/n$','fontsize',fzl,'interpreter','latex')

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
hold on
box on
grid on

if line2==1

gf(1) = plot(itern1,totalarea1-areai1,'color',cl(1,:))
set(gf(1),'linewidth',2);
bf(1) = plot(itern1,totalarea1-areai1,'s','color',cl(1,:))
set(bf(1),'markersize',10,'linewidth',2);

gf(2) = plot(itern2,totalarea2-areai2,'color',cl(2,:))
set(gf(2),'linewidth',2);
bf(2) = plot(itern2,totalarea2-areai2,'d','color',cl(2,:))
set(bf(2),'markersize',10,'linewidth',2);
set(gca,'fontsize',fz);
%set(gca,'xscale','log')

%xlabel('log$(1/n^2)$','fontsize',fzl,'interpreter','latex')
text(9,0.144,'$t_n$ (cutoff time)','fontsize',fzl,'interpreter','latex')
ylabel('$\Delta^3$','fontsize',fzl,'interpreter','latex')
axis([6 14 0.15 0.21])
legend(bf,'Numerical Diffusion','Smoothing Step','location','northeast')
title('$\Delta^3$ versus cutoff time','fontsize',fzl,'interpreter','latex')

else
gf(1) = plot(itern1,totalarea1-areai1,'color',cl(1,:))
set(gf(1),'linewidth',2);
bf(1) = plot(itern1,totalarea1-areai1,'s','color',cl(1,:))
set(bf(1),'markersize',10,'linewidth',2);
set(gca,'fontsize',fz);
text(10.3,0.145,'$t_n$ (cutoff time)','fontsize',fzl,'interpreter','latex')
ylabel('$\Delta^3$','fontsize',fzl,'interpreter','latex')
axis([8 14 0.15 0.20])
%legend(bf,'Numerical Diffusion','Smoothing Step','location','northeast')
title('$\Delta^3$ versus cutoff time','fontsize',fzl,'interpreter','latex')
end







