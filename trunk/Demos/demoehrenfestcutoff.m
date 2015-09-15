%demodhrenfestcutoff.m
%
% Tzu-Chen Liang 11-11-2007

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fz = 14;
fzl = 16;
obj = 0;  %obj=1 slide, obj=0 paper


if obj == 1
 cm = lines(6);
 linestylelist= {'-','-','-','-','-','-'};
 fz=16; fzl=18
else
 cm = zeros(6,3);
 linestylelist= {'-','--','-.','s','d','x'};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


iter = 15000;
dlist = [ 100 200 400 1000 2000];
%dlist= [10000]
dl    = length(dlist);

close all
figure 
hold on
nlist = []

for i=1:dl
 
  n = dlist(i)
  [A,xe] = ehrenfestA(n);
  x0 = zeros(n+1,1);
  xe = binopdf(0:n,n,0.5)';
  x0(1) = 1;
  x = x0;
  nlisti = [];
  for j = 1:iter
     x = A'*x;
     nlisti = [nlisti sum(abs(x-xe))/2];
     %nlist = [nlisti sqrt(sum(xe.*(x./xe).^2))];
  end  
   nlist(i,:) =nlisti;


slist = 1:150:15000;  


end

if obj==0
 for i = 4:dl
    gh(i)= plot(slist, nlist(i,slist));
    set(gh(i),'linewidth',1,'color',cm(i,:),'markersize',12)
 end




for i = 1:dl
   if i<3
     h(i)= plot(nlist(i,:))
   else
     h(i)= plot(slist, nlist(i,slist));
   end
   set(h(i),'linewidth',2,'color',cm(i,:),'linestyle',linestylelist{i},'markersize',6)
end

else
 for i = 1:dl
   h(i)= plot(nlist(i,:))
   set(h(i),'linewidth',2,'color',cm(i,:),'linestyle',linestylelist{i},'markersize',6)
 end
end



wsc = plot([0 4000],[0.5 0.5]);
set(wsc,'linewidth',2,'color',[0 0 0],'linestyle','--');


for i = 1:dl
  w{i} = sprintf('n = %d',dlist(i));
end
legend(h,w,'fontsize',fzl,'interpreter','latex');



grid on
box on
xlabel('$k$ (iteration)','fontsize',fzl,'interpreter','latex')
ylabel('$|\omega^k_n$-$\bar{\omega}|_{TV}$','fontsize',fzl,'interpreter','latex')
title('Random walk on an $n$-dimensional hypercube','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',fz);
axis([0 4000 0 1.4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
     h(i)= plot([0:iter-1]/itern(i),nlist(i,:))
   else
     h(i)= plot(slist/itern(i), nlist(i,slist));
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


legend(h,w,'fontsize',fzl,'interpreter','latex');
axis([0 5 0 1])
set(h,'linewidth',2,'markersize',6);

box on
grid on
xlabel('$k$/$t_n$ (normalized iteration)','fontsize',fzl,'interpreter','latex')
ylabel('$|\omega^k_n$-$\bar{\omega}|_{TV}$','fontsize',fzl,'interpreter','latex')
title('Random walk on an $n$-dimensional hypercube','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',fz);
axis([0 3 0 1.4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intrupthere

iter = 60000;
dlist= [10000]
dl    = length(dlist);


figure 
hold on
nlist = []

for i=1:dl
 
  n = dlist(i)
  [A,xe] = ehrenfestA(n);
  x0 = zeros(n+1,1);
  xe = binopdf(0:n,n,0.5)';
  x0(1) = 1;
  x = x0;
  nlisti = [];
  for j = 1:iter
     x = A'*x;
     nlisti = [nlisti sum(abs(x-xe))/2];
     %nlist = [nlisti sqrt(sum(xe.*(x./xe).^2))];
  end  
   nlist(i,:) =nlisti;


slist = 1:150:15000;  


end

if obj==0
 for i = 4:dl
    gh(i)= plot(slist, nlist(i,slist));
    set(gh(i),'linewidth',1,'color',cm(i,:),'markersize',12)
 end




for i = 1:dl
   if i<3
     h(i)= plot(nlist(i,:))
   else
     h(i)= plot(slist, nlist(i,slist));
   end
   set(h(i),'linewidth',2,'color',cm(i,:),'linestyle',linestylelist{i},'markersize',6)
end

else
 for i = 1:dl
   h(i)= plot(nlist(i,:))
   set(h(i),'linewidth',2,'color',cm(i,:),'linestyle',linestylelist{i},'markersize',6)
 end
end






for i = 1:dl
  w{i} = sprintf('n = %d',dlist(i));
end
legend(h(1),w{1},'fontsize',fzl,'interpreter','latex');



grid on
box on
xlabel('iteration','fontsize',fzl,'interpreter','latex')
ylabel('$|\omega^k_n$-$\bar{\omega}|_{TV}$','fontsize',fzl,'interpreter','latex')
title('Random walk on an $n$-dimensional hypercube','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',fz);
axis([0 60000 0 1.4])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear xi yi areai
figure 
hold on
box on
grid on
xi = 0:1e-4:3;
for i = 1:dl

 yi(i,:) = interp1([0:iter-1]/itern(i),nlist(i,:),xi,'spline');
 areai(i) = trapz(xi,abs(yi(i,:)-0.5));
end

g1 = plot(log(dlist.*log(dlist)/4),3-areai)
set(g1,'linewidth',2);
g2 = plot(log(dlist.*log(dlist)/4),3-areai,'s')
set(g2,'markersize',10,'linewidth',2);
set(gca,'fontsize',fz);



