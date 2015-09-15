function [Xa1,var2,var1,mnorm,snorm,sv,Mm,alln] = testAf(Astruc,Mstruc,ic,nit,showpic,showvar,saveresult)

%
% This function simulates the system by
%  x_k+1 = (A*M^p) * x_k
% with some initial distrubution 
%  
%  ic   = 'cosy','cosx','squx','squy','rand'
%  var1    => 1-norm  
%  var2    => 2 norm
%  mixnorm => mix-norm 
%
% Tzu-Chen Liang     6-16-2006 
%

close all
load lk2000;
Hlist =[0,-0.5,-1];




if isa(Astruc,'struct')
   A  = Astruc.A;
   n  = Astruc.n;
else
   A  = Astruc;
   n  = fix(sqrt(size(A,1)));
   clear Astruc;
   Astruc.mapfunction = 'none';
   Astruc.n = n;
   Astruc.param =[];
end

if ~isempty(Mstruc)
  if isa(Mstruc,'struct')
     M  = Mstruc.M;
     k  = Mstruc.k;
  else
     M  = Mstruc;
  end
else
  k = 1e-3;
end

Xa = zeros(n,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is used to evaluate the mix-norm
%mlist = [n/2:-1:1,0:1:n/2-1];
%nlist = mlist;
%P     = ones(n,1)*mlist.^2+ (ones(n,1)*nlist.^2)';
%P     = fftshift(1./((1+4*pi^2.*P).^0.5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<2
   M=[];
end

if nargin<3
   ic = 'cosx'
end
if nargin<4
   nit = 50;
end
if nargin<5
   showpic = 1;
end
if nargin<6
   showvar = 0;
end
if nargin<7
   saveresult = 0; 
end

    dx = 1/n;
    xs = dx/2 :dx:1-dx/2;
    ys = xs;
    [Xs,Ys] = meshgrid(xs,ys);

if isa(ic,'char')
  if ic == 'cosx'
      Xa = 1*cos(2*pi*Ys);
  elseif ic == 'cosy'
      Xa = 1*cos(2*pi*Xs);
  elseif ic == 'squx'
      Xa(1:fix(n/2),:) = 0.5;
      Xa(fix(n/2)+1:end,:) = -0.5;
  elseif ic == 'squy'
      Xa(:,1:fix(n/2)) = 0.5;
      Xa(:,fix(n/2)+1:end) = -0.5;
  elseif ic == 'rand'
      Xa = rand(n,n);
  elseif ic == 'sing'
      Xa(1:5,1:5) = 1;
  end
else 
    Xa  = ic;
    ic  = 'given';
end

  xa0   = reshape(Xa,n^2,1); 
  xM    = max(xa0);
  xm    = min(xa0);
  xa1   = xa0;
  xave  = ones(n^2,1)*sum(xa1)/n^2;
  xa1   = xa1 - sum(sum(xa1))/n/n;
 

  dy    = 1/n;
  y     = dy/2:dy:1-dy/2;
  z     = dy/2:dy:1-dy/2;
  [Y,Z] = meshgrid(y,z);


  mj      = 1; 
  var2    = [];
  var1    = [];
  mnorm = [];
  snorm   = [];
  dnorm   =[];
  sv = [];

  for i = 1:nit
     
     if i>1         
           if ~isempty(Mstruc)   
             Xa1   = real(ifft2(Xf.*M));
             xa1   = reshape(Xa1,n^2,1);
           end
           xa1 = A*xa1;
     end   
   
   
     % 2-norm
     var2  = [var2 sqrt(var(xa1))];
     % 1-norm
     var1  = [var1 sum(abs(xa1-sum(sum(xa1))/n/n))/n/n];
     %mix-norm
     Xa1 = reshape(xa1,n,n);
     Xf  = fft2(Xa1);
     Xf(1,1) = 0; % make the zero frequency term zero!
     %mnorm = [mixnorm sqrt(sum(sum(P.*(abs(Xf).^2))))/n/n];
     mnorm = [mnorm, mixnorm(Xf,'f',lk)];
     %sobolev norm
     [w] = sobolevnorm(Xf,'f',2,Hlist);
     snorm  = [snorm w'];
     [w] = diffusenorm(Xf,'f',k);
     dnorm  = [dnorm w];   

     % sv
     %sv = [ sv  sum(sum(abs(Xa1-Xinf)))/n/n];
  
     %movie part
     if showpic ==1;          
          %imagesc(Xa1',[xm xM]);
          imagesc(Xa1');
          colormap(gray)
          axis equal 
          axis([1 n 1 n]);
          set(gca,'ytick',[])
          set(gca,'xtick',[])
          box on        
          set(gca,'ydir','normal');
          Mm(mj) = getframe;
          mj = mj+1;
     else
        disp(i)
        Mm=[];
     end      
     

  end
xa1   = xa1 - sum(sum(xa1))/n/n;
Xa1 = reshape(xa1,n,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if showvar==1


   fz = 14;
   fzl = 16;
   obj = 0;  %obj=1 slide, obj=0 paper
   if obj == 1
     cm = jet(6);
     linestylelist= {'-','-','-','-','-','-'};
   else
     cm = zeros(6,3);
     linestylelist= {'-','--','-.','d','x','s'};
   end


  figure
  hold on
  h1(1) = plot(1:nit,var1,'color',cm(1,:),'linestyle',linestylelist{1});
  
  h1(2) = plot(1:nit,mnorm,'color',cm(2,:),'linestyle',linestylelist{2});
  h1(3) = plot(1:nit,snorm(1,:),'color',cm(3,:),'linestyle',linestylelist{3});
  h1(4) = plot(1:nit,snorm(2,:),'color',cm(4,:),'linestyle',linestylelist{4});
  h1(5) = plot(1:nit,snorm(3,:),'color',cm(5,:),'linestyle',linestylelist{5});
  h1(6) = plot(1:nit,dnorm,'color',cm(6,:),'linestyle',linestylelist{6});
  set(h1,'linewidth',2);
    legend('L_1-norm','mix-norm',['H_{',num2str(Hlist(1)),'}-norm'],...
       ['H_{',num2str(Hlist(2)),'}-norm'],['H_{',num2str(Hlist(3)),'}-norm'],'diffuse-norm')
  grid on
  xlabel('iterations','fontsize',fz);
  ylabel('||.||','fontsize',fz);
  set(gca,'fontsize',fz);
  saveas(gcf,'normplot'   ,'eps')


  figure
  hold on
  h2(1) = plot(1:nit,var1/var1(1),'color',cm(1,:),'linestyle',linestylelist{1});
 
  h2(2) = plot(1:nit,mnorm/mnorm(1),'color',cm(2,:),'linestyle',linestylelist{2});
  h2(3) = plot(1:nit,snorm(1,:)/snorm(1,1),'color',cm(3,:),'linestyle',linestylelist{3});
  h2(4) = plot(1:nit,snorm(2,:)/snorm(2,1),'color',cm(4,:),'linestyle',linestylelist{4});
  h2(5) = plot(1:nit,snorm(3,:)/snorm(3,1),'color',cm(5,:),'linestyle',linestylelist{5});
  h2(6) = plot(1:nit,dnorm/dnorm(1),'color',cm(6,:),'linestyle',linestylelist{6});
  set(h2,'linewidth',2);
  legend('$L_1$-norm','mix-norm',['H_{',num2str(Hlist(1)),'}-norm'],...
       ['H_{',num2str(Hlist(2)),'}-norm'],['H_{',num2str(Hlist(3)),'}-norm'],'diffuse-norm')
  set(gca,'YScale','log')
  grid on
  box on
  axis([1 nit 0 1]);
  xlabel('iteration','fontsize',fz)
  ylabel('||.||, normalized','fontsize',fz)
  set(gca,'fontsize',fz);
  saveas(gcf,'lognormplot','eps')
  
end

  alln = [var1/var1(1); 
          mnorm/mnorm(1);
          snorm(1,:)/snorm(1,1);
          snorm(2,:)/snorm(2,1); 
          snorm(3,:)/snorm(3,1);
           dnorm/dnorm(1)];

if saveresult == 1
 paramstr  = num2str(cell2mat(Astruc.param),'%6g');
 while length(strfind(paramstr,'  '))>=1, paramstr = strrep(paramstr,'  ',' ');, end;

 paramstr  = strrep(paramstr,' ','_'); 
 mypath    = ['/home/tzuchen/Desktop/Temp computation result/'];
 fname     = ['testAf_n',num2str(Astruc.n),'_nit',num2str(nit),'_ic_',ic,'_',Astruc.mapfunction,...
           '_',paramstr];
 savestr   = [mypath,fname];

 if ~isempty(Mstruc)
   Mstr    = ['_k',num2str(Mstruc.k),'_dt',num2str(Mstruc.dt),'_period',num2str(Mstruc.period)];
   savestr = [savestr,Mstr];
 end

 save([savestr,'.mat'],'var2','var1','mnorm')


 disp(sprintf('Computation Result is saved'))
end



