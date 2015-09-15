%demo standardmapparam
% need to manually choose 'latex' for the legend
obj=1


if obj == 0
 cm = lines(6);
 linestylelist= {'-','-','-','-','-','-'};
 linestylelist2= {'-','-','-','-','-','-'};
else
 cm = zeros(6,3);
 linestylelist= {'-','--','-.','d','x','s'};
 linestylelist2= {'-','--','-.',':','x','s'};
end

fz  = 14;
fzl = 16;

close all
clear w
figure
hold on
box on
grid on
load mixnorm_40000_501
load mixnorm_40000_503
load mixnorm_40000_505
load mixnorm_40000_507
load mixnorm_40000_509

w(1,:)= mixnorm_40000_501(:,3);
w(2,:)= mixnorm_40000_503(:,3);
w(3,:)= mixnorm_40000_505(:,3);
w(4,:)= mixnorm_40000_507(:,3);
w(5,:)= mixnorm_40000_509(:,3);


for i=1:5
    w(i,:) = w(i,:)/w(i,1);
   
    %h(i)=plot(w(i,:));
    h(i)=plot(w(i,:),'linestyle',linestylelist{i},'linewidth',2,'markersize',7) ;
    set(h(i),'linewidth',2,'color',cm(i,:))    
    
end



set(gca,'fontsize',fz)
axis([0 50 0 1.1])
xlabel('iteration','fontsize',fzl,'interpreter','latex');
ylabel('var$(f)$','fontsize',fzl,'interpreter','latex');
[LEGH,OBJH,OUTH,OUTM] =legend('$\epsilon=0.1$','$\epsilon=0.3$','$\epsilon=0.5$','$\epsilon=0.7$','$\epsilon=0.9$','location','southwest')
set(OBJH(1:5),'interpreter','latex','fontsize',fzl);
title('$f^0=$cos$(2 \pi x_1)$ and $\epsilon=\{0.1, 0.3, 0.5,0.7,0.9\}$','interpreter','latex','fontsize',fzl )

for i=4:5
 thline = plot(w(i,:),'color',[0 0 0],'linewidth',1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear w
figure
hold on
box on
grid on
load mixnorm_40000_501y
load mixnorm_40000_503y
load mixnorm_40000_505y
load mixnorm_40000_507y
load mixnorm_40000_509y

w(1,:)= mixnorm_40000_501y(:,3);
w(2,:)= mixnorm_40000_503y(:,3);
w(3,:)= mixnorm_40000_505y(:,3);
w(4,:)= mixnorm_40000_507y(:,3);
w(5,:)= mixnorm_40000_509y(:,3);


for i=1:5
    w(i,:) = w(i,:)/w(i,1);
   
    %h(i)=plot(w(i,:));
    h(i)=plot(w(i,:),'linestyle',linestylelist{i},'linewidth',2,'markersize',7) ;
    set(h(i),'linewidth',2,'color',cm(i,:))   

end



set(gca,'fontsize',fz)
axis([0 50 0 1.1])
xlabel('iteration','fontsize',fzl,'interpreter','latex');
ylabel('var$(f)$','fontsize',fzl,'interpreter','latex');
[LEGH,OBJH,OUTH,OUTM] =legend('$\epsilon=0.1$','$\epsilon=0.3$','$\epsilon=0.5$','$\epsilon=0.7$','$\epsilon=0.9$','location','southwest')
set(OBJH(1:5),'interpreter','latex','fontsize',fzl);
title('$f^0=$cos$(2 \pi x_2)$ and $\epsilon=\{0.1,0.3,0.5,0.7,0.9\}$','interpreter','latex','fontsize',fzl )

for i=4:5
 thline = plot(w(i,:),'color',[0 0 0],'linewidth',1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear w
figure
hold on
box on
grid on
load mixnorm_40000_8000x
load mixnorm_40000_8000y
load mixnorm_40000_8001x
load mixnorm_40000_8001y

w(1,:) = mixnorm_40000_8000x(:,3);
w(2,:) = mixnorm_40000_8001x(:,3);
w(3,:) = mixnorm_40000_8000y(:,3);
w(4,:) = mixnorm_40000_8001y(:,3);


for i=1:4
    w(i,:) = w(i,:)/w(i,1);
    h(i)=plot(w(i,:),'linestyle',linestylelist2{i},'linewidth',2,'markersize',1) ;
    %h(i)=plot(w(i,:));
    set(h(i),'linewidth',2,'color',cm(i,:));    

end


set(gca,'fontsize',fz);
axis([0 800 0.5 1.05]);
xlabel('iteration','fontsize',fzl,'interpreter','latex');
ylabel('var$(f)$','fontsize',fzl,'interpreter','latex');


[LEGH,OBJH,OUTH,OUTM] =legend('$f^0=$cos$(2\pi x_1),\epsilon=0$','$f^0=$cos$(2\pi x_1),\epsilon=0.1$','$f^0=$cos$(2\pi x_2),\epsilon=0$','$f^0=$cos$(2\pi x_2),\epsilon=0.1$','location','southwest');
%legend(LEGH,'boxoff'\);


set(OBJH(1:4),'interpreter','latex','fontsize',fzl);
title('$\epsilon=\{0,0.1\}$ and $f^0=\{$cos$(2 \pi x_1),$cos$(2 \pi x_2)\}$','interpreter','latex' ,'fontsize',fzl)







