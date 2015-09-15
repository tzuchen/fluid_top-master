
rundata = 0;
if rundata ==1
  [Astruc] = maprefine2(500,[],@standardmap,0.4,1);
  [M]      = fourierdiffuseM(500,1,1e-3,1);
  [Xa1,var2,var1,mixnorm] = testAf(Astruc,[],'cosx',300,0,1,1);
  [Xa1,var2,var1,mixnorm] = testAf(Astruc,[],'cosy',300,0,1,1);
  [Xa1,var2,var1,mixnorm] = testAf(Astruc,[],'squx',300,0,1,1);
  [Xa1,var2,var1,mixnorm] = testAf(Astruc,[],'squy',300,0,1,1);

  [Xa1,var2,var1,mixnorm] = testAf(Astruc,M,'cosx',300,0,1,1);
  [Xa1,var2,var1,mixnorm] = testAf(Astruc,M,'cosy',300,0,1,1);
  [Xa1,var2,var1,mixnorm] = testAf(Astruc,M,'squx',300,0,1,1);
  [Xa1,var2,var1,mixnorm] = testAf(Astruc,M,'squy',300,0,1,1);
end


close all
clear all
mypath   = '/home/tzuchen/Desktop/Temp computation result/'
fname{1} = 'testAf_n500_nit300_ic_cosx_standardmap_0.4_1.mat';
fname{2} = 'testAf_n500_nit300_ic_cosy_standardmap_0.4_1.mat';
fname{3} = 'testAf_n500_nit300_ic_squx_standardmap_0.4_1.mat';
fname{4} = 'testAf_n500_nit300_ic_squy_standardmap_0.4_1.mat';

fname{5} = 'testAf_n500_nit300_ic_cosx_standardmap_0.4_1_k1_dt0.001_period1.mat';
fname{6} = 'testAf_n500_nit300_ic_cosy_standardmap_0.4_1_k1_dt0.001_period1.mat';
fname{7} = 'testAf_n500_nit300_ic_squx_standardmap_0.4_1_k1_dt0.001_period1.mat';
fname{8} = 'testAf_n500_nit300_ic_squy_standardmap_0.4_1_k1_dt0.001_period1.mat';

for i = 1:size(fname,2)
  load([mypath,fname{i}]);
  figure
  hold on
  nit  = length(var1); 
  h(1) = plot(1:nit,var1/var1(1),'r');
  h(2) = plot(1:nit,var2/var2(1),'g');
  h(3) = plot(1:nit,mixnorm/mixnorm(1),'b');
  set(h,'linewidth',2);
  legend('1-norm','2-norm','mix-norm')
  title(fname{i})
  grid on
end


