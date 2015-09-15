% Given A50 and A100
%
% Tzu-Chen 6/1/2006

clear all

mapfunction = @standardmap;
param{1}    = 0.1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n            = 1500;
nit          = 500;
showpic      = 0;


%k            = 1e-3; % real k = k/dt/4/pi^2
dt 	     = 1e-3;
p  	     = 2;
klist       = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70]/1000;
%klist        = [0 5 10 15 20 25 30 35 40]/1000;
%klist        = [0 5 10]/1000;

[A]	     = maprefine2(n,[],mapfunction,param{1:end});
%[M]          = rdwalkM(n,p,k,dt,1);
numofcase    = length(klist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(klist)
  disp(i)
  %[Xa1,var2(i,:)] = testAs(A,M,klist(i),'cosx',nit,showpic);
 
  [M] = fourierdiffuseM(n,klist(i),dt,1);
  [Xa1,var2(i,:), var1(i,:)] = testAf(A,M,p,'cosx',nit,showpic);
end

logvar2 = log(var2);
for i = 1: numofcase
  var2(i,:) = var2(i,:)/var2(i,1);
  vslope(i)  = -(logvar2(i,end) - logvar2(i,end-100+1))/100;
  varn(i,:)  = logvar2(i,:)/vslope(i decay)+[1:nit];    
end


save /home/tzuchen/Desktop/T'emp computation result'/varykdata.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots
  
  sigmaplot       = 1;
  logsigmaplot    = 1;
  asymptplot      = 1;
  normalendplot   = 1;
  normalslopeplot = 1;
  stationary      = 1;
  
  caselist     = klist;
  numofcase    = length(klist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
if sigmaplot ==1
   figure
   hold on

   h = plot(var2(:,1:50)')
   set(h,'linewidth',2)   
   legend(cellstr(num2str(caselist')))
   grid on
   axis([1 50 0 1.05])
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
  h =  plot(varn(1:numofcase,1:end)')
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
    %nylogvar2(i,:) = nlogvar2(i,:)+[1:nit];  
    %nylogvar2(i,:) =nylogvar2(i,:)/nylogvar2(i,end); 
  end

  h= plot(nlogvar2')
  %h= plot(nylogvar2');
  set(h,'linewidth',2)  
  legend(cellstr(num2str(caselist')))
  %axis([1 nit -1.5 0])
  grid on

  title('....')
  xlabel('iterations')
  ylabel('normalized \sigma')


end


