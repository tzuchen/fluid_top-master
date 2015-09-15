

nswitch = 10
[Mn,varx,vary,h0norm,ind] = testltm2(S,[],[],'squy',300,1,200,400,nswitch);
close all
fz = 14;


flist = [1 5 20 40 80 120 160 200]
flist = [1 3 8 12 20 40 60 100]
flist = flist(end:-1:1);

q = -1;
d= 0.01;

for i = 1:length(flist)

f = flist(i);

p  =mod(i,2);
q = q+p;

p ,q
h = axes('Position',[ p/2+0*d q/4+0*d  0.5-d 0.25-d ]);

Mf = Mn(f).cdata;
imagesc(Mf(end:-1:1,:,:));

%imagesc(Mn(f).cdata);
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'TickLength',[0 0]);
grid off
s = sprintf('k=%d',flist(i))
text(0,-10,s,'fontsize',fz)
axis equal
axis tight

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
vary = vary-min(vary);
vary = vary/max(vary);
h= plot(vary)
axis([1 100 0 1])
set(h,'linewidth',2)
grid on
box on 
xlabel('iterations')
ylabel('Total Variation Distance')


