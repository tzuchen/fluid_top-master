function [Astruc] = maprefine2(n,xpi,mapfunction,varargin)
%
% 
%
% Tzu-Chen Liang 6-12-2006
%

 dx       = 1/n; 
 dy       = 1/n;
 x        = dx/2:dx:1-dx/2;
 y        = dy/2:dy:1-dy/2;

 [X,Y]    = meshgrid(x,y);
 % x,y switched!
 yc       = reshape(X,n^2,1);
 xc       = reshape(Y,n^2,1);
 
 param    = varargin;
tic
 [xce,yce]= feval(mapfunction,xc,yc,1,-1,param);

 xe       = 0:dx/2:1;
 ye       = 0:dy/2:1;
 [Ye,Xe]  = meshgrid(xe,ye);

 [Xe,Ye]  = feval(mapfunction,Xe,Ye,0,-1,param);
toc
tic

 xe1      = reshape(Xe(1:2:end-2,1:2:end-2),n^2,1)+1e-5;
 ye1      = reshape(Ye(1:2:end-2,1:2:end-2),n^2,1)+1e-5;

 xe2      = reshape(Xe(2:2:end-1,1:2:end-2),n^2,1);
 ye2      = reshape(Ye(2:2:end-1,1:2:end-2),n^2,1)+1e-5;

 xe3      = reshape(Xe(3:2:end  ,1:2:end-2),n^2,1)-1e-5;
 ye3      = reshape(Ye(3:2:end  ,1:2:end-2),n^2,1)+1e-5;

 xe4      = reshape(Xe(1:2:end-2,2:2:end-1),n^2,1)+1e-5;
 ye4      = reshape(Ye(1:2:end-2,2:2:end-1),n^2,1);

 xe5      = reshape(Xe(2:2:end-1,2:2:end-1),n^2,1);
 ye5      = reshape(Ye(2:2:end-1,2:2:end-1),n^2,1);

 xe6      = reshape(Xe(3:2:end  ,2:2:end  ),n^2,1)-1e-5;
 ye6      = reshape(Ye(3:2:end  ,2:2:end  ),n^2,1);

 xe7      = reshape(Xe(1:2:end-2,3:2:end  ),n^2,1)+1e-5;
 ye7      = reshape(Ye(1:2:end-2,3:2:end  ),n^2,1)-1e-5;

 xe8      = reshape(Xe(2:2:end-1,3:2:end  ),n^2,1);
 ye8      = reshape(Ye(2:2:end-1,3:2:end  ),n^2,1)-1e-5;

 xe9      = reshape(Xe(3:2:end  ,3:2:end  ),n^2,1)-1e-5;
 ye9      = reshape(Ye(3:2:end  ,3:2:end  ),n^2,1)-1e-5;


xp = zeros(9,n^2);
yp = zeros(9,n^2);

 xp(1,:) = mod(ceil(xe1/dx),n);
 xp(2,:) = mod(ceil(xe2/dx),n);
 xp(3,:) = mod(ceil(xe3/dx),n);
 xp(4,:) = mod(ceil(xe4/dx),n);
 xp(5,:) = mod(ceil(xe5/dx),n);
 xp(6,:) = mod(ceil(xe6/dx),n);
 xp(7,:) = mod(ceil(xe7/dx),n);
 xp(8,:) = mod(ceil(xe8/dx),n);
 xp(9,:) = mod(ceil(xe9/dx),n);
 yp(1,:) = mod(ceil(ye1/dy),n);
 yp(2,:) = mod(ceil(ye2/dy),n);
 yp(3,:) = mod(ceil(ye3/dy),n);
 yp(4,:) = mod(ceil(ye4/dy),n);
 yp(5,:) = mod(ceil(ye5/dy),n);
 yp(6,:) = mod(ceil(ye6/dy),n);
 yp(7,:) = mod(ceil(ye7/dy),n);
 yp(8,:) = mod(ceil(ye8/dy),n);
 yp(9,:) = mod(ceil(ye9/dy),n);

yp(find(yp==0))=n;
xp(find(xp==0))=n;



toc
nzind     = zeros(1,9*n^2);
cind      = 1;
ni        = 9;
tic
for i = 1:n^2
       k         = (i-1)*n^2;
       nzind(1,cind:cind+ni-1) = reshape(k+(yp(:,i)-1)*n+xp(:,i),1,ni);
       cind  = cind + ni;
end

toc
 nzind       = nzind(1:cind-1);
 nzind       = unique(sort(nzind));
 nnzA        = length(nzind);
 [jind,iind] = ind2sub([n^2,n^2],nzind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Av = zeros(nnzA,1);
 sumithrow =  zeros(n^2,1);

tic
 for k = 1:nnzA
   i         = iind(k);
   j         = jind(k);
   px        = min(abs(xc(j)-xce(i)),1);
   py        = min(abs(yc(j)-yce(i)),1);   
   Av(k)     = min(abs( dx-min(px,1-px ) )*abs( dy-min(py,1-py))/dx/dy,1);
   sumithrow(i) = sumithrow(i)+Av(k); 
 end
toc
tic  
 for k = 1:nnzA  
   Av(k)     = Av(k)/sumithrow(iind(k));
 end 
toc

 Ap  =sparse(iind,jind,Av,n^2,n^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Astruc.A           = Ap;
Astruc.mapfunction =  func2str(mapfunction);
Astruc.param       = param;
Astruc.n           = n;
Astruc.xpi         = xpi;

 




