%demotentmapcutoffplot
%
% Tzu-Chen Liang 11-11-2007

close all
figure
hold on
box on
grid on
fzl = 24;
if obj == 1
 cm = jet(6);
 linestylelist= {'-','-','-','-','-','-'};
else
 cm = zeros(6,3);
  linestylelist= {'-','-','-','-','-','-'};
end


k=15
nlist  = 3:2:15;

for n=nlist 
  for i =0:k
    if i<=n-1
      mu(i+1) =  1-2^(1+i-n); 
     mu
    else
      mu(i+1) = 0;
    end   
  end

 
   text(n-1.7,mu(n-1)+0.01*n-0.16,sprintf('n=%d',n),'fontsize',18,'interpreter','latex')

   h(i) =  plot(0:k,mu,'linewidth',2)


end
 
text(4,-0.09,'$k$ (iteration)','fontsize',fzl,'interpreter','latex')
ylabel('$|\omega_n^k$-$\bar{\omega}|_{TV}$','fontsize',fzl,'interpreter','latex')
%title('Tent map cutoff','fontsize',fzl,'interpreter','latex')
text(4,1.035,'Tent map cutoff','fontsize',fzl,'interpreter','latex')

set(gca,'fontsize',fz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure 
hold on
box on
grid on

x = 0:0.01:1;
h = plot(x,tentmap(x));
set(h,'linewidth',2)
xlabel('$x$','fontsize',fzl,'interpreter','latex');
ylabel('S$(x)$','fontsize',fzl,'interpreter','latex');
%title('Tent map','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',fz);
%text(0.4,0.4,'Tent map','fontsize',fzl,'interpreter','latex')
text(0.36,1.04,'Tent map','fontsize',fzl,'interpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
hold on
box on
grid on

mu(1) = 1
for i = 1:6
  mu(i+1) = mu(i)/2;
  %stairs([mu(i),1],[1/mu(i),0]);
  h = stairs([0 mu(i)],[1/mu(i),0]);
set(h,'linewidth',2)
 s = ['$\omega_',num2str(i) ,'^0$'];
 text(mu(i),1/mu(i)+1,s,'fontsize',fzl,'interpreter','latex')
end
set(gca,'fontsize',fz);
xlabel('$x$','fontsize',fzl,'interpreter','latex');
ylabel('$\omega(x)$','fontsize',fzl,'interpreter','latex');
text(0.18,36.5, 'Initial distributions $\omega_n^0$','fontsize',fzl,'interpreter','latex')








