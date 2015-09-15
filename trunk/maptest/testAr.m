function [Xa1,var2,var1,mnorm,snorm] = testAr(mapfunction,param,n,Mstruc,ic,nit,showpic,showvar,saveresult)
%
% This function simulates the system by
% the map plus diffution
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




if ~isempty(Mstruc)
  if isa(Mstruc,'struct')
     M  = Mstruc.M;
  else
     M  = Mstruc;
  end
    
end
k   =1e-3;
%n =  size(M,1);
Xa = zeros(n,n);

Astruc.param = param;
Astruc.mapfunction = func2str(mapfunction);
Astruc.n     = n;

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
  end
else 
    Xa  = ic;
end
  Xf    = fft2(Xa);
  xa0   = reshape(Xa,n^2,1); 
  xM    = max(xa0);
  xm    = min(xa0);
  xa1   = xa0;
  xave  = ones(n^2,1)*sum(xa1)/n^2;

  dy    = 1/n;
  y     = dy/2:dy:1-dy/2;
  z     = dy/2:dy:1-dy/2;
  [Y,Z] = meshgrid(y,z);

  yl     = -dy/2:dy:1+dy/2;
  zl     = -dy/2:dy:1+dy/2;
  [Yl,Zl] = meshgrid(yl,zl);

  [Yi,Zi]=feval(mapfunction,Y,Z,1,-1,param);

  mj      = 1; 
  var2    = [];
  var1    = [];
  mnorm = [];
  snorm   = [];
  dnorm   =[];

  for i = 1:nit
     
     if i>1         

 
        if isempty(Mstruc)
           
              if ic == 'cosx'
                 Xa1 = 1*cos(2*pi*Yi);
              elseif ic == 'cosy'
                 Xa1 = 1*cos(2*pi*Zi);
              elseif ic == 'squx'
                 Xa1(find(Yi>0.5))  = 0.5;
                 Xa1(find(Yi<=0.5)) = -0.5;
              elseif ic == 'squy'
                 Xa1(find(Zi>0.5)) = 0.5;
                 Xa1(find(Zi<=0.5)) = -0.5;
              end
           if ~isempty(Mstruc)  
             Xf    = fft2(Xa1); 
             Xa1   = real(ifft2(Xf.*M));
             xa1   = reshape(Xa1,n^2,1);
           end

         else
            Xal = [Xa1(end,end) Xa1(end,:) Xa1(end,1);
                   Xa1(:,end)   Xa1        Xa1(:,1)  ;
                   Xa1(1,end)   Xa1(1,:)   Xa1(1,1)  ];

            Xa1 = interp2(Yl,Zl,Xal,Yi,Zi);
            if ~isempty(Mstruc)  
               Xf    = fft2(Xa1); 
               Xa1   = real(ifft2(Xf.*M));
               xa1   = reshape(Xa1,n^2,1);
            end
         end

        [Yi,Zi]=feval(mapfunction,Yi,Zi,1,-1,param);

           xa1 = reshape(Xa1,n^2,1); 
     end   
       
     % 2-norm
     var2  = [var2 sqrt(var(xa1))];
     % 1-norm
     var1  = [var1 sum(abs(xa1-xave))/n/n];
     %mix-norm
     Xa1 = reshape(xa1,n,n);
   
     %mnorm = [mixnorm sqrt(sum(sum(P.*(abs(Xf).^2))))/n/n];
     mnorm = [mnorm, mixnorm(Xf,'f',lk)];
     %sobolev norm
     [w] = sobolevnorm(Xf,'f',2,[0,-0.5,-1]);
     snorm  = [snorm w'];
     [w] = diffusenorm(Xf,'f',k);
     dnorm  = [dnorm w];    
     %movie part
     if showpic ==1;          
          imagesc(Xa1,[xm xM]);
          axis equal   
          axis([1 n 1 n]);
          box on
          set(gca,'ytick',[])
          set(gca,'xtick',[])
          set(gca,'ydir','normal');
          Mm(mj) = getframe;
          mj = mj+1;
     else
        disp(i)
     end      
     

  end
 
Xa1 = reshape(xa1,n,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mnorm

if showvar==1
  figure
  hold on
  h(1) = plot(1:nit,var1/var1(1),'r');
  h(2) = plot(1:nit,var2/1,'g');
  h(3) = plot(1:nit,mnorm,'b');
  h4 = plot(1:nit,snorm);
  h(4) = plot(1:nit,dnorm,'y');
  set(h,'linewidth',2);
  legend('1-norm','2-norm','mix-norm','sobolevnorm')
  grid on
end

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

