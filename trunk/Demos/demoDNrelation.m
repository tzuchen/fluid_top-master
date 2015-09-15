%demoDNrelation

mapfunction = @standardmap;
param{1}    = 0.2;
ic          = 'cosx';
iter        = 150;

Nlist = [50 100 200 500 1000];
var2  = [];

for i = 1:length(Nlist)

  n = Nlist(i);
  [Astruc] = maprefine2(n,[],mapfunction,param{:});
  Dstar(i) = Destimate(mapfunction,param,n);
  [Mdiff]  = fourierdiffuseM(n,Dstar(1)-Dstar(i),1,1);
  [Xa{i},var2i,var1,mnorm,w,sv] = testAf(Astruc,Mdiff,ic,iter,0,1,0);
  var2(i,:) = var2i;



end

save /home/tzuchen/Desktop/T'emp computation result'/demoDNrelationdata
load /home/tzuchen/Desktop/T'emp computation result'/demoDNrelationdata

close all

h = semilogy(var2');
set(h,'linewidth',2)
grid on
box on


figure
for i = 1:5
subplot(3,2,i)
imagesc(Xa{i}')
axis equal

end


