% Demo example2plotfull
% Tzu-Chen Liang 11-07-2007
close all
%clear all

%load example2plotfulldata
%load Sstroock
%load S50full
%load S100full
%load S200full
%load S400full
%load S800full

fzl= 16;
fz = 14;
obj = 1;
nit    = 30;
U = 1.2 %flow velocity
l = 0.01   %characteristic length

nlist = [20 50 100 200 400];

tlist = [0:nit-1]*0.015*8;
lx     = 0.015*8

Dlist = [ 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7]*1e-4


%Dlist = (0.01./nlist).^2 *U/lx;
Pelist = U*l./Dlist; 


if obj == 1
 cm = lines(6);
 linestylelist= {'-','-','-','-','-','-'};
else
 cm = zeros(6,3);
 linestylelist= {'-','--','-.','d','x','s'};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear MS mydata


%[xa1,MS{1},mydata{1}] = testA2(S800,Dlist(1)/1e-4*lx,800,1600);
%[xa1,MS{2},mydata{2}] = testA2(S800,Dlist(2)/1e-4*lx,800,1600);
%[xa1,MS{3},mydata{3}] = testA2(S800,Dlist(3)/1e-4*lx,800,1600);
%[xa1,MS{4},mydata{4}] = testA2(S800,Dlist(4)/1e-4*lx,800,1600);
%[xa1,MS{5},mydata{5}] = testA2(S800,Dlist(5)/1e-4*lx,800,1600);
%[xa1,MS{6},mydata{6}] = testA2(S800,Dlist(6)/1e-4*lx,800,1600);

lxs = 0.01*14;
%[xa1,Mk{1},mdata{1}] = testA2(Sstroock,Dlist(1)/1e-4*lxs,800,1600);
%[xa1,Mk{2},mdata{2}] = testA2(Sstroock,Dlist(2)/1e-4*lxs,800,1600);
%[xa1,Mk{3},mdata{3}] = testA2(Sstroock,Dlist(3)/1e-4*lxs,800,1600);
%[xa1,Mk{4},mdata{4}] = testA2(Sstroock,Dlist(4)/1e-4*lxs,800,1600);
%[xa1,Mk{5},mdata{5}] = testA2(Sstroock,Dlist(5)/1e-4*lxs,800,1600);
%[xa1,Mk{6},mdata{6}] = testA2(Sstroock,Dlist(6)/1e-4*lxs,800,1600);


lx3d = 0.04
U3d = 0.18;
Pelist3d = U3d*l./Dlist;
t3dlist = [0:nit-1]*0.02;
%[xa1,M3{1},m3d{1}] = testA2(S3d,Dlist(1)/1e-4*lx3d,800,1600);
%[xa1,M3{2},m3d{2}] = testA2(S3d,Dlist(2)/1e-4*lx3d,800,1600);
%[xa1,M3{3},m3d{3}] = testA2(S3d,Dlist(3)/1e-4*lx3d,800,1600);
%[xa1,M3{4},m3d{4}] = testA2(S3d,Dlist(4)/1e-4*lx3d,800,1600);
%[xa1,M3{5},m3d{5}] = testA2(S3d,Dlist(5)/1e-4*lx3d,800,1600);
%[xa1,M3{6},m3d{6}] = testA2(S3d,Dlist(6)/1e-4*lx3d,800,1600);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close 
figure 
hold on
box on
grid on
for i=1:6

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

for i = 1:6
  casen{i} = sprintf('Pe = %2.2e',Pelist(i));
end

legend(h,casen)

for i = 4:6
  gh(i) =plot(tlist,mydata{i},'linewidth',1,'color',cm(i,:));
end


plot([0 3],[ 0.05 0.05],'linewidth',2,'linestyle','--')
text(0.2,0.066,'$x_{90}$','fontsize',fzl,'interpreter','latex')
set(gca,'fontsize',fz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stroockPe = [2e3 2e4 2e5 9e5];
stroocky90 = [0.66 0.88 1.32 1.7];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on


for i = 1:6
   y90(i)   = min(find(mydata{i}<0.05))*lx;
   y90s(i)   = min(find(mdata{i}<0.05))*lxs;
   y903d(i) = min(find(m3d{i}<0.05))*lx3d;
end

g(1) = semilogx(Pelist,y90,'color',cm(1,:));
g2(1) = plot(Pelist,y90,'s','markersize',8,'linewidth',2,'color',cm(1,:));

g(2) = plot(Pelist3d,y903d,'linewidth',2,'color',cm(4,:));
g2(2) = plot(Pelist3d,y903d,'o','markersize',8,'linewidth',2,'color',cm(4,:));

g(3)= plot(stroockPe,stroocky90,'linewidth',2,'color',cm(2,:));
g2(3) = plot(stroockPe,stroocky90,'d','markersize',8,'linewidth',2,'color',cm(2,:));


g(4) = plot(Pelist,y90s,'linewidth',2,'color',cm(3,:));
g2(4) = plot(Pelist,y90s,'x','markersize',8,'linewidth',2,'color',cm(3,:));

set(gca,'xscale','log')
set(g,'linewidth',2);
set(gca,'fontsize',fz)
xlabel('Pe','fontsize',fzl,'interpreter','latex')
ylabel('$x_{90}$(cm)','fontsize',fzl,'interpreter','latex')
title('Mixing length versus Pe','fontsize',fzl,'interpreter','latex')
%axis([10^(2.5) 10^6.1 0.8 2])
box on
grid on

legend(g2,'optimal herringbone structure','optimal 3-D structure', 'Stroock''s experiments','simulation of Stroock''s  structure','location','northwest')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on


g(1) = plot(log(Pelist),y90,'color',cm(1,:));
g2(1) = plot(log(Pelist),y90,'s','markersize',8,'linewidth',2,'color',cm(1,:));

g(2) = plot(log(Pelist3d),y903d,'linewidth',2,'color',cm(4,:));
g2(2) = plot(log(Pelist3d),y903d,'o','markersize',8,'linewidth',2,'color',cm(4,:));

g(3) = plot(log(stroockPe),stroocky90,'linewidth',2,'color',cm(2,:));
g2(3) = plot(log(stroockPe),stroocky90,'d','markersize',8,'linewidth',2,'color','r','color',cm(2,:));

g(4) = plot(log(Pelist),y90s,'linewidth',2,'color',cm(3,:));
g2(4) = plot(log(Pelist),y90s,'x','markersize',8,'linewidth',2,'color',cm(3,:));


%set(gca,'xscale','log')
set(g,'linewidth',2);
set(gca,'fontsize',fz)
xlabel('log(Pe)','fontsize',fzl,'interpreter','latex')
ylabel('$x_{90}$(cm)','fontsize',fzl,'interpreter','latex')
title('Mixing length versus log(Pe)','fontsize',fzl,'interpreter','latex')
%axis([10^(2.5) 10^6.1 0.8 2])
box on
grid on
legend(g2,'optimal herringbone structure','optimal 3-D structure', 'Stroock''s experiments','simulation of Stroock''s  structure','location','northwest')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
     Mn = MS{3};
     flist = [ 2 3 4 5 10]

  else
     Mn = MS{6};
     flist = [ 2 3 4 5 10]

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
s = sprintf('x=%2.2fcm',(flist(i)-1)*0.015*8)
text(0,-10,s,'fontsize',fz)
axis equal
axis tight

end



end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pol =1;
if pol==1
figure

q = -1;
d= 0.01;
p = 0;

Mn(1) = M3{4}(16);
Mn(2) = M3{4}(4);
Mn(3) = MS{3}(6);
Mn(4) = Mk{3}(6);



Mn(1).cdata =Mn(1).cdata(:,end:-1:1,:);
Mn(2).cdata =Mn(2).cdata(:,end:-1:1,:);

distancelist = [ 15*lx3d, 3*lx3d 5*lx, 5*lxs,]
cyclelist      = [15 3 5 5];
timelist     = [15*lx3d/0.18, 3*lx3d/0.18, 5*lx/1.2 ,5*lxs/1.2  ]
caselist     =[3 4 1 2]
f = 1
typename={'optimal 3-D structure','optimal 3-D structure','optimal herringbone structure','Stroock''s herringbone design'}

for q =0:1 
  for p= 0:1 
  
 


  h = axes('Position',[ p/2+0*d q/2+0*d  0.5-d 1/2-d ]);

  Mf = Mn(f).cdata;
  imagesc(Mf);

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'TickLength',[0 0]);
grid off
s = sprintf('(%s)%s, x=%2.2fcm, %dth cycle, t= %1.1fsec',char(96+caselist(f)),typename{f},distancelist(f),cyclelist(f),timelist(f));
text(0,-10,s,'fontsize',fz)
axis equal
axis tight
f=f+1

end



end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
q = -1;
d= 0.01;
p = 0; 

flist = [16 6 5 4 3 2];
flist = [5 4 3 2];
for q = 0:3

   h = axes('Position',[ p/3+0*d 0.95*(q/4+0*d)  1/3-d 0.95*(1/4-d) ]);
  
  myf=M3{3}(flist(q+1)).cdata;
 % [fx,fy,fz]=size(myf);
 %  myf = myf(1:ceil(0.75)*fx,:,:);
  imagesc(myf(end:-1:1,end:-1:1,:));

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'TickLength',[0 0]);
grid off
s = sprintf('x=%2.2fcm',(flist(q+1)-1)*lx3d);
text(0,-15,s,'fontsize',fz)
axis equal
axis tight
f=f+1
 if q==3
  s = sprintf('Optimal 3-D Mixer');
  %text(0,245,s,'fontsize',fz)
  text(0,-45,s,'fontsize',fz)
 end

end
p=1
for q = 0:3

    h = axes('Position',[ p/3+0*d 0.95*(q/4+0*d)  1/3-d 0.95*(1/4-d) ]);
  
  myf=MS{3}(flist(q+1)).cdata;
 % [fx,fy,fz]=size(myf);
 %  myf = myf(1:ceil(0.75)*fx,:,:);
  imagesc(myf);

  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'TickLength',[0 0]);
  grid off
s = sprintf('x=%2.2fcm',(flist(q+1)-1)*lx);
text(0,-15,s,'fontsize',fz)
  axis equal
  axis tight
  f=f+1
 if q==3
  s = sprintf('Optimal Harringbone Mixer');
  %text(0,245,s,'fontsize',fz)
  text(0,-45,s,'fontsize',fz)
 end

end

p=2
for q = 0:3


   h = axes('Position',[ p/3+0*d 0.95*(q/4+0*d)  1/3-d 0.95*(1/4-d) ]);
  
  myf=Mk{3}(flist(q+1)).cdata;
 % [fx,fy,fz]=size(myf);
 %  myf = myf(1:ceil(0.75)*fx,:,:);
  imagesc(myf);

  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'TickLength',[0 0]);
  grid off
s = sprintf('x=%2.2fcm',(flist(q+1)-1)*lxs);
text(0,-15,s,'fontsize',fz)
  axis equal
  axis tight
  f=f+1
 
 if q==3
  s = sprintf('Stroock''s Mixer');
  %text(0,245,s,'fontsize',fz)
  text(0,-45,s,'fontsize',fz)
 end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

flist = [4 3 2 1 8 7 6 5 12 11 10 9];
f= 1
for p = 0:2
for q = 0:3

   h = axes('Position',[ p/3+0*d 0.95*(q/4+0*d)  1/3-d 0.95*(1/4-d) ]);
  
  myf=MS{4}(flist(f)).cdata;
 % [fx,fy,fz]=size(myf);
 %  myf = myf(1:ceil(0.75)*fx,:,:);
  imagesc(myf);

  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'TickLength',[0 0]);
  grid off
s = sprintf('x=%2.2fcm',(flist(f)-1)*lx);
text(0,-15,s,'fontsize',fz)
  axis equal
  axis tight
  f=f+1
 
 if and(q==3,p==0)
  s = sprintf('Optimal Herringbone Mixer');
  %text(0,245,s,'fontsize',18)
  text(0,-55,s,'fontsize',18)
 end

end

end


% need a structure plot!


