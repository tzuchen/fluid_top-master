function [Xa1,varx,vary,h0norm,ind] = testltm(A1,A2,M,ic,nit,showmovie)
% 
% This function simulates the system by
%  x_k+1 = (A*M^p) * x_k
% with some initial distrubution 
% 
% Tzu-Chen Liang     4-7-2006 
%

close all
A = A1;




n = sqrt(size(A,1));
Xa = zeros(n,n);

if nargin<3
   M = [];
end
if nargin<6
   showmovie = 0;
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
    ic  = 'given';
end

  xa0   = reshape(Xa,n^2,1); 
  xM    = max(xa0);
  xm    = min(xa0);
  xa1   = xa0;
  xave  = ones(n^2,1)*sum(xa1)/n^2;
  xa1   = xa1 - sum(sum(xa1))/n/n;

xc1 = 0.4;
yc1 = 0.4;
r1  = 0.28;
xc2 = 0.6;
yc2 = 0.6;
r2  = 0.28;

[th,r]= cart2pol(Xs-xc1,Ys-yc1);
ind1 = find(r<r1);  
[th,r]= cart2pol(Xs-xc2,Ys-yc2);
ind2 = find(r<r2);  
ind  = union(ind1,ind2);



  xM = max(xa0);
  xm = min(xa0);

  xa1 = xa0;
  Xa1 = Xa;

  xave = ones(n*n,1)*sum(xa1)/n^2;

  dy = 1/n;
  y  = dy/2:dy:1-dy/2;
  z  = dy/2:dy:1-dy/2;
  [Y,Z] = meshgrid(y,z);

  j = 1;
   
  
   varx = [];
   vary = [];
   h0norm = [];
 
    for i = 1:nit
%disp(i)
        if i>1                      
           if mod(i,30)<15
                A = A1;
           else
                A = A2;
           end
               
           xa1 = A*xa1;
        end      
     
          Xa1 = reshape(xa1,n,n);
        if ~isempty(M)
    
          Xf  = fft2(Xa1);
          Xa1   = real(ifft2(Xf.*M.M)); 
          xa1 = reshape(Xa1,n^2,1);
        
        end

        varx  = [varx sqrt(var(Xa1(ind)))];
        vary  = [vary norm(xa1-0.5)];

        %movie part
        if showmovie ==1
          imagesc(Xa1,[xm xM]);
          colormap(gray)
          axis equal    
          axis tight     
          set(gca,'ydir','normal');
          Mm(j) = getframe;
          j = j+1;
        end
    end


 
x=[];
Xa1 = reshape(xa1,n,n);
