function [Xa1,Mm,mydata] = testA(A,D,sa1,sa2,nit)
% 
% This function simulates the system by
%  x_k+1 = A * x_k
% with some initial distrubution 
% 
% Tzu-Chen Liang     3-24-2006 
%
if nargin<5
nit=30
end

close all
x0=[];
OPTS.disp = 0;

% Only give A
if nargin == 2
    sa = sqrt(size(A,1));
    Xa = zeros(sa,sa);
    dx = 1/sa;
    xs = dx/2 :dx:1-dx/2;
    ys = xs;
    [Xs,Ys] = meshgrid(xs,ys);
    Xa = 1*cos(2*pi*Xs);
    %Xa = 1*cos(2*pi*Ys);
    Xa(1:fix(sa/2),:) = 0;
    Xa(fix(sa/2):end,:) = 1;
    
    

    Xa(:,1:fix(sa/2)) = 0;
    Xa(:,fix(sa/2):end) = 1;
    %Xa = 1*cos(2*pi*Xs);
    xa0 = reshape(Xa,(sa)*(sa),1);
    sa1 = sa;
    sa2 = sa;  
else
    dy = 1/sa1;
    ys = dy/2 :dy:1-dy/2;
    dz = 1/sa2; 
    zs = dz/2 :dz:1-dz/2;
    [Ys,Zs] = meshgrid(ys,zs);
    Xa = zeros(sa1,sa2);

    Xa(1:fix(sa1/2),:) = 0;
    Xa(fix(sa1/2):end,:) = 1;
    
    Xa(:,1:fix(sa2/2)) = 1;
    Xa(:,fix(sa2/2):end) = 0;
    
  
    xa0 = reshape(Xa,(sa1)*(sa2),1); 



% Initial distribution is given
 %   if prod(size(x0))>1;     
 %      xa0 = x0;
 %   end




end
  xM = max(xa0);
  xm = min(xa0);

  xa1 = xa0;

  dy = 1/sa1;
  y  = dy/2:dy:1-dy/2;
  dz = 1/sa2;
  z  = dz/2:dz:1-dz/2;
  [Y,Z] = meshgrid(y,z);

  j = 1;
  xave = ones(sa1*sa2,1)*sum(xa1)/(sa1*sa2);  
  
  sig   = [];
  sig1  = [];

 
  [Mstruc] = fourierdiffuseM(max(sa1,sa2),D,1,1);
   M = Mstruc.M([1:sa1/2,end-sa1/2+1:end],:) ;
 

    for i = 1:nit




        if i>1  
           xa1 = A*xa1;
        end      
 
        sig1 = [sig1 1/2*sum(abs(xa1-xave))/sum(xa0) ];
        sig = [sig sqrt(var(xa1)) ];

        if ismember(i,[1:1:1400])        
           Xa1 = reshape(xa1,sa1,sa2);

      

          %for qq= 1:6
	  %  Xa1 = diffu(Xa1);
	  %end
        if i>1
          Xf = fft2(Xa1);
          Xa1   = real(ifft2(Xf.*M));
          xa1   = reshape(Xa1,sa1*sa2,1);
        end
         
        Xt = Xa1([0.15*sa1:0.85*sa1],[0.15*sa2:0.85*sa2]);
        %Xt = Xa1([0.2*sa1:0.8*sa1],[0.2*sa2:0.8*sa2]);
        [sa1r,sa2r] = size(Xt);  
        mydata(i) = sqrt(var(reshape(Xt,sa1r*sa2r,1)));  


          imagesc(Xa1(end:-1:1,end:-1:1),[xm xM]);
          xa1 = reshape(Xa1,sa1*sa2,1);
         
           % Xaw = [Xa1 Xa1;Xa1 Xa1];
           %imagesc(Xaw,[xm xM]);
           
  
           
           %Xa1 = diffu(Xa1);
           %Xa1 = diffu(Xa1);
           %Xa1 = diffu(Xa1);
           %xa1 = reshape(Xa1,sa1*sa2,1);
           
   

           set(gca,'XTick',[],'XTickLabel',[])
           set(gca,'YTick',[],'YTickLabel',[])
           
           axis equal 
           %axis([1 sa2 1 sa1]) 
           axis tight
           colormap gray 
          % pause        
           Mm(j) = getframe;
           j = j+1;
        end

    end
 %save Mm  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   [V,D]       = eigs(A,6,'LM',OPTS);
%   [dr,I]      = sort(diag(abs(D)));
   
  
        
% d1  = dr(end);
% d2  = dr(end-1);
% d3  = dr(end-2);
% iter = 1:nit;
% d1list = sig(1);
% d2list = sig(1);
% d3list = sig(1);

% for i = 2:nit
%   d1list = [d1list d1*d1list(end)];
%   d2list = [d2list d2*d2list(end)];
%   d3list = [d3list d3*d3list(end)];
% end

% d2rlist =[sig(end)];
% for i = 2:nit
%   d2rlist = [d2rlist(1)/d2 d2rlist];
% end


 %mydata = [sig1;sig;d1list;d2list;d2rlist];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xf = diffu(X)

C0 = 1/8;
C1 = 1/4;
C2 = 3/16;
 
[sa1,sa2] = size(X);

X0 = X;
Xt = zeros(sa1,sa2);

Xf = C0*X;

Xt(1:end-1,:) = X(2:end,:);
Xf = Xf+C1*Xt;
Xt = Xt*0;

Xt(2:end,:) = X(1:end-1,:);
Xf = Xf+C1*Xt;
Xt = Xt*0;

Xt(1:end-2,:) = X(3:end,:);
Xf = Xf+C2*Xt;
Xt = Xt*0;

Xt(3:end,:) = X(1:end-2,:);
Xf = Xf+C2*Xt;
Xt = Xt*0;

X=Xf; 

Xf = C0*X;

Xt(:,1:end-1) = X(:,2:end);
Xf = Xf+C1*Xt;
Xt = Xt*0;

Xt(:,2:end) = X(:,1:end-1);
Xf = Xf+C1*Xt;
Xt = Xt*0;

Xt(:,1:end-2) = X(:,3:end);
Xf = Xf+C2*Xt;
Xt = Xt*0;

Xt(:,3:end) = X(:,1:end-2);
Xf = Xf+C2*Xt;


Xf(1:2,:) = X0(1:2,:);
Xf(end-1:end,:) = X0(end-1:end,:);
Xf(:,1:2) = X0(:,1:2);
Xf(:,end-1:end) = X0(:,end-1:end);
Xf(1:2,:) = X0(1:2,:);
Xf(:,end-1:end) = X0(:,end-1:end);

































  

