function [Mm,varx,vary,h0norm,ind] = testltm(A1,A2,M,ic,nit,showmovie,na1,na2,nswitch)
% 
% This function simulates the system by
%  x_k+1 = (A*M^p) * x_k
% with some initial distrubution 
% 
% Tzu-Chen Liang     4-7-2006 
%

close all
A = A1;

if nargin<6
 na1 = sqrt(size(A,1));
 na2 = na1;
end


if isempty(A2)
   p = 1:na1*na2; 
   ps = reshape(p,na1,na2);
   %pr = ps(:,end:-1:1);
   pr = ps(end:-1:1,end:-1:1);
   
   A2 = A(pr,pr);
    %A2 = A(:,pr);

end



Xa = zeros(na1,na2);

if nargin<3
   M = [];
end
if nargin<6
   showmovie = 0;
end

    dx = 1/na1;
    xs = dx/2 :dx:1-dx/2;
    dy = 1/na2;
    ys = dy/2 :dy:1-dy/2;
   
    [Xs,Ys] = meshgrid(xs,ys);

if isa(ic,'char')
  if ic == 'cosx'
      Xa = 1*cos(2*pi*Ys);
  elseif ic == 'cosy'
      Xa = 1*cos(2*pi*Xs);
  elseif ic == 'squx'
      Xa(1:fix(na1/2),:) = 0.5;
      Xa(fix(na1/2)+1:end,:) = -0.5;
  elseif ic == 'squy'
      Xa(:,1:fix(na2/2)) = 0.5;
      Xa(:,fix(na2/2)+1:end) = -0.5;
  elseif ic == 'rand'
      Xa = rand(na1,na2);
  end
else 
    Xa  = ic;
    ic  = 'given';
end

  xa0   = reshape(Xa,na1*na2,1); 
  xM    = max(xa0);
  xm    = min(xa0);
  xa1   = xa0;
  xave  = ones(na1*na2,1)*sum(xa1)/na1/na2;
  xa1   = xa1 - sum(sum(xa1))/na1/na2;

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

  xave = ones(na1*na2,1)*sum(xa1)/na1/na2;

  dy = 1/na1;
  y  = dy/2:dy:1-dy/2;
  dz = 1/na2;
  z  = dz/2:dz:1-dz/2;
  [Y,Z] = meshgrid(y,z);

  j = 1;
   
  
   varx = [];
   vary = [];
   h0norm = [];
 
    for i = 1:nit
%disp(i)
        if i>1                      
           if mod(i,nswitch)<nswitch/2
                A = A1;
           else
                A = A2;
           end
               
           xa1 = A*xa1;
        end      
     
          Xa1 = reshape(xa1,na1,na2);
        if ~isempty(M)
    
          Xf  = fft2(Xa1);
          Xa1   = real(ifft2(Xf.*M.M)); 
          xa1 = reshape(Xa1,na1*na2,1);
        
        end
 
       Xc = Xa1(0.2*na1:0.8*na1,0.2*na2:0.8*na2);
       varx  = [varx sqrt(var(Xc(1:end)))];
       
       % varx  = [varx sqrt(var(Xa1(1:end)))];
        %varx  = [varx sqrt(var(Xa1(ind)))];
        vary  = [vary norm(xa1-0.5)];

        %movie part
        if showmovie ==1
          imagesc(Xa1,[xm xM]);
          colormap(gray)
          axis equal    
          axis tight     
          set(gca,'ydir','normal');
          set(gca,'xtick',[]);
          set(gca,'ytick',[]);
          grid off 
          box on
          Mm(j) = getframe;
          j = j+1;
        end
    end


 
x=[];
Xa1 = reshape(xa1,na1,na2);
