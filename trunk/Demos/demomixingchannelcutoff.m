close all

obj=0;
if obj == 1
 cm = lines(6);
 linestylelist= {'-','-','-','-','-','-'};
else
 cm = zeros(6,3);
 linestylelist= {'-','--','-.','d','x','s'};
end

dl = 5
Pelist = Pelist(1:dl);

for i = 1:dl
  casen{i} = sprintf('Pe = %2.2e',Pelist(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close 
figure 
hold on
box on
grid on
for i=1:dl

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
title('Chage of Peclet Number','fontsize',fzl,'interpreter','latex')

%for i = 1:dl
%  casen{i} = sprintf('Pe = %2.2e',Pelist(i));
%end

legend(h,casen)

for i = 4:dl
  gh(i) =plot(tlist,mydata{i},'linewidth',1,'color',cm(i,:));
end


plot([0 3],[ 0.05 0.05],'linewidth',2,'linestyle','--')
text(0.2,0.066,'$x_{90}$','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',fz)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




W = cell2mat(mydata');
nlist = W;

iter =30
noval = 0.25

for i = 1:dl
  for j = 1:iter-1
    if nlist(i,j+1)>=nlist(i,j);
      nlist(i,j+1) = nlist(i,j)-1e-9; 
    end
  end
end

clear itern
for i = 1:dl
  itern(i)   = interp1(nlist(i,:),0:iter-1,noval);
end

clist = 'bgrcmyk'

figure
hold on


 if obj==0

  for i = 4:dl
    gh(i)= plot([0:iter-1]/itern(i), nlist(i,:));
    set(gh(i),'linewidth',1,'color',cm(i,:))
 end

 for i = 1:dl
   if i<3
     h(i)= plot([0:iter-1]/itern(i),nlist(i,:));
   else
    % h(i)= plot(slist/itern(i), nlist(i,slist)); 
     h(i)= plot([0:iter-1]/itern(i),nlist(i,:));
   end
   set(h(i),'linewidth',2,'color',cm(i,:),'linestyle',linestylelist{i},'markersize',6)
 end

 else

  for i = 1:dl
     h(i) = plot([0:iter-1]/itern(i),nlist(i,:));
     set(h(i),'color',cm(i,:));  
  
   end

end






plot([0 1],[1 1],'linewidth', 2, 'color',[0 0 0],'linestyle','--')
plot([1 1],[1 0],'linewidth', 2, 'color',[0 0 0],'linestyle','--')
plot([1 3],[0 0],'linewidth', 2, 'color',[0 0 0],'linestyle','--')


legend(h,casen,'fontsize',fzl,'interpreter','latex');
axis([0 5 0 0.5])
set(h,'linewidth',2,'markersize',6);

box on
grid on
xlabel('$k$/$t_n$ (normalized iteration)','fontsize',fzl,'interpreter','latex')
ylabel('var$(f^k)$','fontsize',fzl,'interpreter','latex')
title('Chaotic Mixing Channel','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',fz);
%axis([0 3 0 1.4])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear areai;
xi = 0:1e-4:3;
for i = 1:dl

 yi(i,:) = interp1([0:iter-1]/itern(i),nlist(i,:),xi,'spline');
 areai(i) = trapz(xi,abs(yi(i,:)-noval));
end



figure 
hold on
box on
grid on 
 itern1 = itern;
 nlist1 = nlist;
 noval1 = noval; 
 totalarea1 = abs(nlist(end,end)-noval)*3; 
 areai1 = areai;
 %Dlist1 = dlist;

Dlist1= Pelist;
itern1 = itern;
cf(1) = plot(Dlist1,itern1)
set(cf(1),'linewidth',2,'color',cm(1,:));
df(1) = plot(Dlist1,itern1,'s')
set(df(1),'markersize',10,'linewidth',2,'color',cm(1,:));
set(gca,'xscale','log')
set(gca,'fontsize',fz)
xlabel('Pe','fontsize',fzl,'interpreter','latex')
%text(10^8.5,5.25,'$n$','fontsize',fzl,'interpreter','latex')
ylabel('$t_n$ (cutoff time)','fontsize',fzl,'interpreter','latex')
axis([10^(3.5) 10^(8.5)  1 9])
%legend(df,'Numerical Diffusion','Smoothing Step')
title('Cutoff time versus $n$','fontsize',fzl,'interpreter','latex')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
hold on
box on
grid on


  gf(1) = plot(Dlist1,totalarea1-areai1,'color',cm(1,:))
  set(gf(1),'linewidth',2);
  bf(1) = plot(Dlist1,totalarea1-areai1,'s','color',cm(1,:))
  set(bf(1),'markersize',10,'linewidth',2);
  set(gca,'fontsize',fz);
  set(gca,'xscale','log')
  text(10^-8.5,0.144,'$1/n$','fontsize',fzl,'interpreter','latex')
  ylabel('$\Delta$','fontsize',fzl,'interpreter','latex')
  %axis([10^(3.5) 10^(8.5) 0.12 0.24])
  %legend(bf,'Numerical Diffusion','Smoothing Step','location','northwest')
  title('$\Delta$ versus $1/n$','fontsize',fzl,'interpreter','latex')
  xlabel('Pe','fontsize',fzl,'interpreter','latex')









