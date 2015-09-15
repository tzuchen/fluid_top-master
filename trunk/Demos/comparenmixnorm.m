clear all

%mapfunction = @twistmap;
%param{1}    = 0.5;
%param{2}    = 0.5;
%param{3}    = 0.1;

mapfunction = @standardmap;
param{1}    = 0.1;

k            = 1e-4; % real k = k/dt/4/pi^2
dt 	     = 1e-3;
p            = 2;

nlist       = [50 100 150 200 400 600 800 1000 ];

nit         = 500;
showpic     = 0;




numofcase   = length(nlist);
for i = 1:numofcase
   n      = nlist(i) 
   [A{i}] = maprefine2(n,[],mapfunction,param{1:end});
  
    [M]          = fourierdiffuseM(n,k,dt,1);
    [Xa1,var2t,var1t,mixnormt] = testAf(A{i},M,1,'cosx',nit,showpic);

   var2(i,:) = var2t;
   var1(i,:) = var1t;  
   mixnorm(i,:) = mixnormt;

end

var2n = var2/var2(1,1);
var1n = var1/var1(1,1);
mixnormn = mixnorm/mixnorm(1,1);


save /home/tzuchen/Desktop/T'emp computation result'/comparenmixnormdatawithdiff.mat
