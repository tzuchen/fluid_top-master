
clear all



mapfunction  = @arnoldscatmap;
param{1}     = 1e-3;

%nlist       = [50 100 150 ]%;200 400 600 800 1000 1200 1400 1600];
nlist       = [50 100 150 200 400 600 800 1000 1200 1400 1600];

nit         = 20;
showpic     = 0;
showvar     = 0;



numofcase   = length(nlist);
for i = 1:numofcase
   n      = nlist(i) 
   [A{i}] = maprefine2(n,[],mapfunction,param{1:end});
 
    [Xa1,var2t,var1t] = testAf(A{i},[],'cosx',nit,showpic,showvar,0);
                        
   var2(i,:) = var2t;
   var1(i,:) = var1t;  

end

save /home/tzuchen/Desktop/T'emp computation result'/comparen2data.mat
cutoffnormalize(var2,0.3)
close all
