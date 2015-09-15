% Demo example2plot
% Tzu-Chen Liang 11-07-2007
close all
%clear all
%load example2Sdata
%load example2SdataS0
%load example2SdataS1

fzl= 16;
fz = 14;
obj = 0;
nit    = 151;
U = 1.2 %flow velocity
l = 0.01   %characteristic length

nlist = [20 50 100 200 400 800];
ncycle =[6 8 12 16 20 26 ];
tlist = [0:nit-1]*0.015;

lx     = 0.015

Dlist = (0.01./nlist).^2 *U/lx;
Pelist = U*l./Dlist; 


if obj == 1
 cm = jet(6);
 linestylelist= {'-','-','-','-','-','-'};
else
 cm = zeros(6,3);
 linestylelist= {'-','--','-.','d','x','s'};
end

[Mm{1},varx{1},vary,h0norm,ind]=testltm2(S20 ,[],[],'squy',nit,1,20,40,8);
[Mm{2},varx{2},vary,h0norm,ind]=testltm2(S50 ,[],[],'squy',nit,1,50,100,8);
[Mm{3},varx{3},vary,h0norm,ind]=testltm2(S100,[],[],'squy',nit,1,100,200,8);
[Mm{4},varx{4},vary,h0norm,ind]=testltm2(S200,[],[],'squy',nit,1,200,400,8);
[Mm{5},varx{5},vary,h0norm,ind]=testltm2(S400,[],[],'squy',nit,1,400,800,8);
[Mm{6},varx{6},vary,h0norm,ind]=testltm2(S800,[],[],'squy',nit,1,800,1600,8);


[Mb{1},varxb{1},vary,h0norm,ind]=testltm2(S200,[],[],'squy',nit,1,200,400,6);
[Mb{2},varxb{2},vary,h0norm,ind]=testltm2(S200,[],[],'squy',nit,1,200,400,8);
[Mb{3},varxb{3},vary,h0norm,ind]=testltm2(S200,[],[],'squy',nit,1,200,400,12);
[Mb{4},varxb{4},vary,h0norm,ind]=testltm2(S200,[],[],'squy',nit,1,200,400,16);
[Mb{5},varxb{5},vary,h0norm,ind]=testltm2(S200,[],[],'squy',nit,1,200,400,20);
[Mb{6},varxb{6},vary,h0norm,ind]=testltm2(S200,[],[],'squy',nit,1,200,400,26);

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flist = [8 16 24 32 64]-1;
flist = flist(end:-1:1);
Mn = Mm{6};
q = -1;
d= 0.01;
p = 0;

for p =0:1 
  q = -1;
  if p == 0
     Mn = Mm{3};
  else
     Mn = Mm{6};
  end
  for i = 1:length(flist)
  f = flist(i);
  q=q+1


  h = axes('Position',[ p/2+0*d q/5+0*d  0.5-d 1/5-d ]);

  Mf = Mn(f).cdata;
%imagesc(Mf(end:-1:1,:,:));
imagesc(Mf);
%imagesc(Mn(f).cdata);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'TickLength',[0 0]);
grid off
s = sprintf('x=%2.2fcm',flist(i)*0.015)
text(0,-10,s,'fontsize',fz)
axis equal
axis tight

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%flist = [8 16 24 32 40 48 56 64]-1;
%flist = flist(end:-1:1);
%Mn = Mm{3};
%q = -1;
%d= 0.01;

%for i = 1:length(flist)

%f = flist(i);

%p  =mod(i,2);
%q = q+p;


%h = axes('Position',[ p/2+0*d q/4+0*d  0.5-d 0.25-d ]);

%Mf = Mn(f).cdata;
%imagesc(Mf(end:-1:1,:,:));
%imagesc(Mf);
%imagesc(Mn(f).cdata);
%set(gca,'xtick',[]);
%set(gca,'ytick',[]);
%set(gca,'TickLength',[0 0]);
%grid off
%s = sprintf('k=%d',flist(i))
%text(0,-10,s,'fontsize',fz)
%axis equal
%axis tight

%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
box on
grid on

for i = 1:6

  h(i) = plot(tlist, varx{i});
  set(h(i),'color',cm(i,:) );
  set(h(i),'linestyle',linestylelist{i},'linewidth',2);

end

axis([0 2 0 0.55 ])
xlabel('$x$(cm)','fontsize',fzl,'interpreter','latex')
ylabel('Standard Deviation','fontsize',fzl,'interpreter','latex')
title('Chage of Peclet Number','fontsize',fzl,'interpreter','latex')

for i = 1:6
  casen{i} = sprintf('Pe = %2.2d',Pelist(i));

end

legend(casen)

plot([0 2],[ 0.05 0.05],'linewidth',2,'linestyle','--')
text(0.2,0.066,'$x90$','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',fz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
box on
grid on
for i = 1:6

  h(i) = plot(tlist, varxb{i});
  set(h(i),'color',cm(i,:) );
  set(h(i),'linestyle',linestylelist{i},'linewidth',2);

end

axis([0 2 0 0.55 ])
xlabel('$x$(cm)','fontsize',fzl,'interpreter','latex')
ylabel('Standard Deviation','fontsize',fzl,'interpreter','latex')
title('Change of cycles','fontsize',fzl,'interpreter','latex')

for i = 1:6
  casen{i} = sprintf('%d-cycle',ncycle(i));
end

legend(casen,'location','southwest')



set(gca,'fontsize',fz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on

for i = 1:6
   y90(i) = min(find(varx{i}<0.05))*lx;
end

%g = semilogx(Pelist,y90);
%g2 = plot(Pelist,y90,'s','markersize',8,'linewidth',2);
%set(gca,'xscale','log')

g = plot(log(Pelist),y90);
g2 = plot(log(Pelist),y90,'s','markersize',8,'linewidth',2);



set(g,'linewidth',2);
set(gca,'fontsize',fz)
xlabel('Pe','fontsize',fzl,'interpreter','latex')
ylabel('$x90$(cm)','fontsize',fzl,'interpreter','latex')
title('Mixing distance versus Pe','fontsize',fzl,'interpreter','latex')
%axis([10^(2.5) 10^6.1 0.8 2])
box on
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %[xa1,MS0,mydata] = testA(S1,100,200);
 %[xa1,MS0,mydata] = testA(S100,100,200);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pol =1;
if pol==1
figure

q = -1;
d= 0.01;
p = 0;

for p =0:1 
  q = -1;
  if p == 0
     [Mn]=testltm2(S100,[],[],'squy',nit,1,100,200,9);
     %Mn = Mb{3};
     flist = [1 2 3 4 5]*9-1;

  else
     Mn = MS0;
     flist = [2 3 4 5 6]

  end
  flist = flist(end:-1:1);
  for i = 1:length(flist)
  f = flist(i);
  q=q+1


  h = axes('Position',[ p/2+0*d q/5+0*d  0.5-d 1/5-d ]);

  Mf = Mn(f).cdata;
%imagesc(Mf(end:-1:1,:,:));
imagesc(Mf);
%imagesc(Mn(f).cdata);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'TickLength',[0 0]);
grid off
s = sprintf('x=%2.2fcm',flist(i)*0.015)
text(0,-10,s,'fontsize',fz)
axis equal
axis tight

end



end
end


% need a structure plot!






