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

nlist       = [50 100 150 200 400 600 800 1000 1200 1400 1600];
%nlist       = [50 100 150 200];% 400 600 800 1000 1200 1400 1600];
nit         = 500;
showpic     = 0;
showvar     = 0;



numofcase   = length(nlist);
for i = 1:numofcase
   n      = nlist(i) 
   [A{i}] = maprefine2(n,[],mapfunction,param{1:end});
   %[M]          = rdwalkM(n,p,k,dt,1); 
   %[Xa{i},var2t,var1t] = testAs(A{i},[],0,'cosx',nit,showpic);
  
    [M]          = fourierdiffuseM(n,k,dt,1);
    [Xa1,var2t,var1t] = testAf(A{i},[],'cosx',nit,showpic,showvar,0);
                        
   var2(i,:) = var2t;
   var1(i,:) = var1t;  

end

save /home/tzuchen/Desktop/T'emp computation result'/comparendata.mat

close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
  
  sigmaplot       = 1;
  logsigmaplot    = 1;
  asymptplot      = 1;
  normalendplot   = 1;
  normalslopeplot = 1;
  
  caselist = nlist;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if sigmaplot ==1
   figure
   hold on

   h = plot(var2(:,1:nit)')
   set(h,'linewidth',2)   
   %legend(cellstr(num2str(caselist')))
   grid on
   axis([1 nit 0 1.05])
   title('\sigma vs iteration')
   xlabel('iterations')
   ylabel('\sigma')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if logsigmaplot ==1
   figure
   hold on

   for i = 1: numofcase
      var2(i,:) = var2(i,:)/var2(i,1);
   end

   logvar2 = log(var2);

   for i = 1: numofcase
      vslope(i)  = -(logvar2(i,end) - logvar2(i,end-100+1))/100;
      varn(i,:)  = logvar2(i,:)./(vslope(i)*[1:nit]);    
   end

   h =  semilogy(var2(:,1:50)')
   legend(cellstr(num2str(caselist')))
   set(gca,'YScale','log')
   grid on
   set(h,'linewidth',2)   
   title('log(\sigma) vs iteration')
   xlabel('iterations')
   ylabel('log(\sigma)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if asymptplot == 1;

  figure
  hold on
  h =  plot(varn(:,1:end)')
  legend(cellstr(num2str(caselist')))
  grid on
  set(h,'linewidth',2)   
  title('....')
  xlabel('iterations')
  ylabel('log(\sigma)/log(\sigma^{\infty})') 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if normalendplot == 1;
  
  figure
  hold on

  pt2fit = nit-10:nit;
  for i = 1: numofcase
    fitpoly = polyfit(pt2fit,logvar2(i,pt2fit),1);
    fitdata(i,:) = polyval(fitpoly,1:nit);
   
  end

  h= plot(-(logvar2./fitdata)')
  set(h,'linewidth',2)  
  legend(cellstr(num2str(caselist')))
  axis([1 nit -1.5 0])
  grid on

  title('....')
  xlabel('iterations')
  ylabel('log(\sigma)/fittet straight line')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if normalslopeplot ==1 

 
  figure
  hold on

  pt2fit = nit-100:nit;
  for i = 1: numofcase
    fitpoly = polyfit(pt2fit,logvar2(i,pt2fit),1);
    slope(i)  = fitpoly(1);
    %fitdata(i,:) = polyval(fitpoly,1:nit);
    nlogvar2(i,:) = -logvar2(i,:)/slope(i);   
    nylogvar2(i,:) = -nlogvar2(i,:)/nlogvar2(i,end); 

  end

  h= plot(nlogvar2');
  %h= plot(nylogvar2');
  set(h,'linewidth',2)  
  legend(cellstr(num2str(caselist')))
  %axis([1 nit -1.5 0])
  grid on

  title('....')
  xlabel('iterations')
  ylabel('normalized \sigma')


end








