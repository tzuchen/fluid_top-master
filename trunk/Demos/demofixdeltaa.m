%  demo of the total variaiton change 
% 
%
%  Tzu-Chen Liang 10-1-2007 
close all
figure 
hold on

% adjust font size
fz = 14;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


obj = 0;  %obj=1 slide, obj=0 paper


if obj == 1
 cm = jet(6);
 linestylelist= {'-','-','-','-','-','-'};
 cm = zeros(6,3);
 linestylelist= {'-','--','-.',':','d','s'};
 fz = 16;
 fzl = 18;
else
 cm = zeros(6,3);
 linestylelist= {'-','--','-.',':','d','s'};
 fz = 14;
 fzl = 16;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w{1} = '-';
w{2} = '--';
w{3} = '-.';
w{4} = ':';

 deltalist = [ 0.6 0.7 0.8 0.9];

 for j = 1:length(deltalist)
     for i = 1:300
 
 	p = 300-i;
 	delta = deltalist(j);
	ep = delta-0.5;

	ks=  floor(p* (log(2)+log(0.5-ep) )/(log(0.5-ep)-log(0.5+ep)));

	Z = p-ks;
 	W = ks+1;
 	tv(i)  = betainc(0.5,Z,W)- betainc(1-delta,Z,W); 
 
      end
 
      h(j) = plot([-299:0],tv) ;
      set(h(j),'linewidth',2,'linestyle', linestylelist{j}); 
      grid on
      box on
      
      xlabel('-$p$','fontsize',fzl,'interpreter','latex');
      ylabel('$|\omega^0$-$\bar{\omega}|_{TV}$','fontsize',fzl,'interpreter','latex');    
      title('$p$ versus total variation distance','fontsize',fzl,'interpreter','latex');  

 end
s = {'$\delta^M=0.6$','$\delta^M=0.7$','$\delta^M=0.8$','$\delta^M=0.9$'};
legend(h,s,'location','southwest','fontsize',fzl,'interpreter','latex');
set(gca,'fontsize',fz);  
axis auto
set(gca,'PlotBoxAspectRatio',[1.2 1 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
hold on
clear h
r1 = 0.5;
r2 = 0.7;



plist = [500,200, 50];
for i = 1:length(plist)

   p = plist(i);
   x =  0:p;
   y1 = pdf('bino',x,p,r1);
   h(i)= plot(x/p,y1,'r');
   set(h(i),'linewidth',2,'linestyle',linestylelist{i});

end

 s = {'$p=500$','$p=200$','$p=50$'};
legend(h,s,'fontsize',fzl,'interpreter','latex','location','northeast');
xlabel('$i/p$','fontsize',fzl,'interpreter','latex');
ylabel('$\omega(i/p)$','fontsize',fzl,'interpreter','latex');

for i = 1:length(plist)

   p = plist(i);
   x =  0:p;
   y2 = pdf('bino',x,p,r2);
   g(i)= plot(x/p,y2,'b');
   set(g(i),'linewidth',2,'linestyle',linestylelist{i});

end
set(gca,'fontsize',fz);
title('The binomial distributions of different $p$','fontsize',fzl,'interpreter','latex')

grid on
box on
set(gca,'PlotBoxAspectRatio',[1.2 1 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
hold on
clear h
p = 500;
r = 1-1/p;
eps = 0.5
k =[ 250,450, 1000]
%r = [0.9,0.7,0.55];
d = 0.5+eps*r.^k;


x =  0:p;
for i = 1:length(k)

y1 = pdf('bino',x,p,d(i));
h(i) = plot(x,y1);
set(h(i),'linewidth',2,'linestyle',w{i});
end

s = {'$k=250$','$k=500$','$k=1000$'};

%s = {'$\delta^M=0.9$','$\delta^M=0.7$','$\delta^M=0.45$'};
legend(h,s,'fontsize',fzl,'interpreter','latex');

y1 = pdf('bino',x,p,0.5);
hi = plot(x,y1);
set(hi,'linewidth',2,'color','r','linestyle','-');

set(gca,'fontsize',fz);
xlabel('$i$','fontsize',fzl,'interpreter','latex');
ylabel('$\omega(i)$','fontsize',fzl,'interpreter','latex');

grid on
box on
set(gca,'PlotBoxAspectRatio',[1.2 1 1])
title('The binomial distributions of different $k$','fontsize',fzl,'interpreter','latex')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
figure 
hold on
clear h
rlist = [0.998];
slist = 1:100:3000

for i = 1:length(rlist)

   r= rlist(i);
   delta0 = 0.9999;
   delta = delta0;
   ep = delta-0.5;
   for k = 1:3000
  
      p = 1/(1-r);
      p = floor(p);
  
      ks=  floor(p* (log(2)+log(0.5-ep) )/(log(0.5-ep)-log(0.5+ep)));
      if(isnan(ks)), ks = 0; end

      Z = p-ks;
      W = ks+1;
      tvu(k)  = betainc(0.5,Z,W)- betainc(1-delta,Z,W); 

      epm = ep*r^p;
      deltam= epm+0.5;
      ks=  floor(p* (log(2)+log(0.5-epm) )/(log(0.5-epm)-log(0.5+epm)));
      if(isnan(ks)), ks = 0; end
      Z = p-ks;
      W = ks+1;
      tvl(k)  = betainc(0.5,Z,W)- betainc(1-deltam,Z,W); 


      %tvue(k) = erf(sqrt(p)*ep);
      %tvle(k) = erf(sqrt(p)*epm);

      tvue(k) = erf(sqrt(p)*ep/sqrt(2));
      tvle(k) = erf(sqrt(p)*epm/sqrt(2));


      %yu = pdf('bino',0:p,p,ep+0.5);       
      %yl = pdf('bino',0:p,p,epm+0.5);       
      %ym = pdf('bino',0:p,p,0.5);       
   
      %tvut(k) = sum(abs(yu-ym))/2;
      %tvlt(k) = sum(abs(yl-ym))/2;


      ep = ep*r;
      delta = ep+0.5;
   end

       h(1) = plot(tvu,'linestyle',linestylelist{1}) ;
       set(h(1),'linewidth',2); 
       h(2) = plot(tvl,'r','linestyle',linestylelist{2}) ;
       set(h(2),'linewidth',2) ;

       h(3) = plot(slist,tvue(slist),'linestyle',linestylelist{5}) ;
       set(h(3),'linewidth',2); 
       h(4) = plot(slist,tvle(slist),'r','linestyle',linestylelist{6}) ;
       set(h(4),'linewidth',2) ;
 
      % h(5) = plot(tvut,'linestyle',linestylelist{5}) ;
      % set(h(5),'linewidth',2); 
      % h(6) = plot(tvlt,'r','linestyle',linestylelist{6}) ;
      % set(h(6),'linewidth',2) ;


       grid on
       box on
       xlabel('iteration','fontsize',fzl,'interpreter','latex');
       ylabel('$|\omega_n^k$-$\bar{\omega}|_{TV}$','fontsize',fzl,'interpreter','latex');  
       title('Upper bound and lower bound','fontsize',fzl,'interpreter','latex'); 
       s = {'Upper bound $(\delta^M, p)$','Lower bound $(\delta^m, p)$','approx. Ub','approx. Lb'}

       legend(h,s,'fontsize',15,'interpreter','latex' );
       set(gca,'fontsize',fz);
end

set(gca,'PlotBoxAspectRatio',[1.2 1 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
hold on
clear h
box on
rlist = [0.998];
slist = 1:100:3000


 h(1) = plot(tvu,'linestyle',linestylelist{1}) ;
       set(h(1),'linewidth',2); 
 plot([250 500 1000],tvu([250 500 1000]),'s','markersize',8,'linewidth',2);
 axis([0 3000 0 1.05])
 

grid on
xlabel('$k$(iteration)','fontsize',fzl,'interpreter','latex');
ylabel('$|\omega_n^k$-$\bar{\omega}|_{TV}$','fontsize',fzl,'interpreter','latex'); 
title('Total variation versus iteration','fontsize',fzl,'interpreter','latex');
set(gca,'fontsize',fz);
 
