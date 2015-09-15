% Demo add diffusion
% 
% Tzu-Chen Liang 7-28-2006 

%n           = 600;
cutratio    = 0.01;
iter        = 60;
ic          = 'cosx'
mapfunction = @standardmap;
param{1}    = 0.2
%mapfunction = @arnoldscatmap;
%param{1}    = 1e-3;


Dlist = [1e-5,1e-4,1e-3,1e-2];

nlist = Nestimate(mapfunction,param,Dlist);
n =200;%nlist(1);
Dstar = Destimate(mapfunction,param,n);


[P,x0,statelist,X0] = freqmap(n,cutratio,iter,ic,mapfunction,param);
x0 = x0/norm(x0);
normx = [];


for k  = 1: length(Dlist)
  Pd = freqmapaddD(P,n,statelist,Dstar,Dlist(k));
  x = x0;
  for i= 1:iter
     normx(k, i) = norm(x);
     x = Pd'*x; 
  end

end

 n = nlist(1);
 var2 = [];
 for i = 1:length(nlist);

   [Mdiff]  = fourierdiffuseM(n,Dlist(i)-Dstar,1,1);   
   [Astruc] = maprefine2(n,[],mapfunction,param{:});
   [Xa1,var2(i,:),var1,mnorm,snorm,sv] = testAf(Astruc,Mdiff,ic,iter,0,0,0);
   var2(i,:) = var2(i,:)/var2(i,1);
 end


save /home/tzuchen/Desktop/T'emp computation result'/demoDdata
load /home/tzuchen/Desktop/T'emp computation result'/demoDdata







close all
figure
hold on
h = plot(normx')
set(h,'linewidth',2);
h = plot(var2','--')
set(h,'linewidth',2);
axis([1 60 0 1 ])

mylegend = char(reshape(double(sprintf('D=%0.0e, ',Dlist')),9,length(Dlist))')

legend(mylegend);
hold off
grid on
box on

figure
hold on
h = semilogy(normx')
set(h,'linewidth',2);
h = semilogy(var2','--')
set(h,'linewidth',2);
set(gca,'YScale','log')
axis([1 25 1e-3 1 ])
mylegend = char(reshape(double(sprintf('D=%0.0e, ',Dlist')),9,length(Dlist))')
legend(mylegend);
grid on
box on




