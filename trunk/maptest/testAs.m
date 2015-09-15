function [Xa1,var2,var1,var3,var4] = testAs(A,M,p,ic,nit,showpic)

%
% This function simulates the system by
%  x_k+1 = (A*M^p) * x_k
% with some initial distrubution 
%  
%  ic   = 'cosx','squx','squy'
%  var1 => 1-norm  
%  var2 => 2 norm
% 
% Tzu-Chen Liang     6-16-2006 
%

close all
sa = sqrt(size(A,1));
Xa = zeros(sa,sa);


if nargin<2
   M=[];
end
if nargin<3
   p = 0;
end
if nargin<4
   ic = 'cosx'
end
if nargin<5
   nit = 50;
end
if nargin<6
   showpic = 1;
end

    dx = 1/sa;
    xs = dx/2 :dx:1-dx/2;
    ys = xs;
    [Xs,Ys] = meshgrid(xs,ys);

if isa(ic,'char')
  if ic == 'cosx'
      Xa = 1*cos(2*pi*Ys);
  elseif ic == 'cosy'
      Xa = 1*cos(2*pi*Xs);
  elseif ic == 'squx'
      Xa(1:fix(sa/2),:) = 1;
      Xa(fix(sa/2)+1:end,:) = 0;
  elseif ic == 'squy'
      Xa(:,1:fix(sa/2)) = 0;
      Xa(:,fix(sa/2)+1:end) = 1;
  elseif ic == 'sing'
      Xa(1:100,1:100) = 1;
  elseif ic == 'strp'
      %Xa(1:125,:) = ones(125,1)*cos(2*pi*xs);
      Xa(1:125,:) = ones(125,sa)/125;
      %Xa=Xa';
  end
else 
    Xa  = ic;
end

  xa0   = reshape(Xa,(sa)*(sa),1); 
  xM    = max(xa0);
  xm    = min(xa0);
  xa1   = xa0;
  xave  = ones(sa*sa,1)*sum(xa1)/sa^2;

  dy    = 1/sa;
  y     = dy/2:dy:1-dy/2;
  z     = dy/2:dy:1-dy/2;
  [Y,Z] = meshgrid(y,z);

  mj    = 1; 
  var2  = [];
  var1  = [];
  var3  = [];
  var4  = []; 

  for i = 1:nit

     if i>1         
        for j = 1:p
           xa1 = M*xa1;
        end
           %xt = sum(reshape(xa1,sa,sa));
           %xa1 = reshape((ones(sa,1)*xt),sa^2,1);
                 %xa1= reshape(reshape(xa1,sa,sa)',sa*sa,1);

           xa1 = A*xa1;
     end   
       
     % 2-norm
     var2  = [var2 sqrt(var(xa1))];
     % infinity-norm
     var1  = [var1 max(max(abs(xa1)))];

     %movie part
     if showpic ==1;
          Xa1 = reshape(xa1,sa,sa);
          imagesc(Xa1,[xm xM]); 
         %  imagesc(Xa1)
    var3=[var3  sum(abs(sum(Xa1')-sum(sum(Xa1'))/sa))];
    var4=[var4  sum(abs(sum(Xa1)-sum(sum(Xa1))/sa))];

%  plot(sum(Xa1')-sum(sum(Xa1'))/sa)

          box on          
          axis equal  
          axis tight       
          set(gca,'ydir','normal');
          set(gca,'XTick',[],'XTickLabel',[])
          set(gca,'YTick',[],'YTickLabel',[])
          Mm(mj) = getframe;
          mj = mj+1;
     else
          disp(i);
     end      

  end
 
Xa1 = reshape(xa1,sa,sa);

