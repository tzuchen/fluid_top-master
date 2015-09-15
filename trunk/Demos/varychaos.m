clear all

mapfunction = @standardmap;
n           = 100;
nit         = 100;
ic          = 'cosx';
k           = 1e-3;
dt          = 1;
period      = 1;
showpic     = 0;
showvar     = 0;
saveresult  = 0;

chaoslist   = [0  0.2  0.4  0.6  0.8  1];

numofcase   = length(chaoslist);


[M]         = fourierdiffuseM(n,k,dt,period);  % no use here
xpi         = ones(n^2,1)/n^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numofcase
   [A{i}] = maprefine2(n,xpi,mapfunction,chaoslist(i));
   [Xa{i},var2t,var1t,mixnormt] = testAf(A{i}',[],ic,nit,showpic,showvar,saveresult);
   var2(i,:)    = var2t;
   var1(i,:)    = var1t; 
   mixnorm(i,:) = mixnormt;
end

save /home/tzuchen/Desktop/T'emp computation result'/varychaosdata.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:numofcase
   [A{i}] = maprefine2(n,xpi,mapfunction,chaoslist(i));
   [Xa{i},var2t,var1t,mixnormt] = testAf(A{i}',M,ic,nit,showpic,showvar,saveresult);
   var2(i,:)    = var2t;
   var1(i,:)    = var1t; 
   mixnorm(i,:) = mixnormt;
end

save /home/tzuchen/Desktop/T'emp computation result'/varychaosdataM.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot part
load /home/tzuchen/Desktop/T'emp computation result'/varychaosdata.mat
close all
var1n    = var1/var1(1,1);
var2n    = var2/var2(1,1);
mixnormn = mixnorm/mixnorm(1,1);
clear h
figure
hold on
h1 = plot(var1n','linewidth',2);
%h2 = plot(var2n','linewidth',2);
h3 = plot(mixnormn','linewidth',2);
legend(num2str(chaoslist'))
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
load /home/tzuchen/Desktop/T'emp computation result'/varychaosdataM.mat
close all
var1n    = var1/var1(1,1);
var2n    = var2/var2(1,1);
mixnormn = mixnorm/mixnorm(1,1);
clear h
figure
hold on
h1 = plot(var1n','linewidth',2);
%h2 = plot(var2n','linewidth',2);
h3 = plot(mixnormn','linewidth',2);
grid on








