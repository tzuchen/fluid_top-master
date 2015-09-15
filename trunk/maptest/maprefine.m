function [Ap,Bp] = maprefine(n,xpi,mapfunction,varargin)

%
% Tzu-Chen Liang 6-12-2006
%

 %n = fix(sqrt(size(A,1)));



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

 xe       = 0:dx:1;
 ye       = 0:dy:1;
 [Ye,Xe]  = meshgrid(xe,ye);


 [Xe,Ye]  = feval(mapfunction,Xe,Ye,0,-1,param);

toc
tic
 xe1      = reshape(Xe(1:end-1,1:end-1),n^2,1);
 ye1      = reshape(Ye(1:end-1,1:end-1),n^2,1);
 xe2      = reshape(Xe(2:end  ,1:end-1),n^2,1);
 ye2      = reshape(Ye(2:end  ,1:end-1),n^2,1);
 xe3      = reshape(Xe(1:end-1,2:end  ),n^2,1);
 ye3      = reshape(Ye(1:end-1,2:end  ),n^2,1);
 xe4      = reshape(Xe(2:end  ,2:end  ),n^2,1);
 ye4      = reshape(Ye(2:end  ,2:end  ),n^2,1);


 xl       =  ceil(min(min(xe1/dx,xe3/dx),min(xe2/dx,xe4/dx))+1*1e-15);
 xr       =  ceil(max(max(xe2/dx,xe4/dx),max(xe1/dx,xe3/dx))-1*1e-15);
 yd       =  ceil(min(min(ye3/dy,ye4/dy),min(ye1/dy,ye2/dy))+1*1e-15);
 yu       =  ceil(max(max(ye1/dy,ye2/dy),max(ye3/dy,ye4/dy))-1*1e-15);


toc
nzind     = zeros(1,10*n^2);
cind      = 1;
tic
for i = 1:n^2
       k         = (i-1)*n^2;
       xlist     = mod(xl(i):xr(i),n);
       ylist     = mod(yd(i):yu(i),n);
       xlist(find(xlist==0))=n;
       ylist(find(ylist==0))=n;


       [Xi,Yi]   = meshgrid(xlist,ylist);
       ni        = prod(size(Xi));
       %nzind     = [ nzind ,reshape(k+(Yi-1)*n+Xi,1,ni)];

        nzind(1,cind:cind+ni-1) = reshape(k+(Yi-1)*n+Xi,1,ni);
        cind  = cind + ni;

end


toc
nzind = nzind(1:cind-1);


 nzind   = unique(sort(nzind));
 nnzA    = length(nzind);
 [jind,iind] = ind2sub([n^2,n^2],nzind);

 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Av = zeros(nnzA,1);
 sumithrow =  zeros(n^2,1);
 sumjthcol =  zeros(n^2,1);
tic
 for k = 1:nnzA
   i = iind(k);
   j = jind(k);
   px = min(abs(xc(j)-xce(i)),1);
   py = min(abs(yc(j)-yce(i)),1);   

   Av(k) = min(abs( dx-min(px,1-px ) )*abs( dy-min(py,1-py))/dx/dy,1);
   sumithrow(i) = sumithrow(i)+Av(k); 
   sumjthcol(j) = sumjthcol(j)+Av(k);
 end
toc
tic  
Bv = Av; 
 for k = 1:nnzA  
   Av(k) = Av(k)/sumithrow(iind(k));
   Bv(k) = Bv(k)/sumjthcol(jind(k));
 end 
toc

 Ap  =sparse(iind,jind,Av,n^2,n^2);
 Bp  =sparse(iind,jind,Bv,n^2,n^2);
 




