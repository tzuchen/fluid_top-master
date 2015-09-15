function [A] = map2markov(n,xpi,mapfunction,varargin)
% Given a map, this function returns the Markov matrix
% A with the same invariant distribution xpi
% 
% mapfunction has to be a function handle
% see standardmap.m for example 
% It has to be periodical and has period 1. 
%
% Tzu-Chen Liang 3-24-2006
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
 [xce,yce]= feval(mapfunction,xc,yc,1,-1,param);

 xe       = 0:dx:1;
 ye       = 0:dy:1;
 [Ye,Xe]  = meshgrid(xe,ye);

 [Xe,Ye]  = feval(mapfunction,Xe,Ye,0,-1,param);

 xe1      = reshape(Xe(1:end-1,1:end-1),n^2,1);
 ye1      = reshape(Ye(1:end-1,1:end-1),n^2,1);
 xe2      = reshape(Xe(2:end  ,1:end-1),n^2,1);
 ye2      = reshape(Ye(2:end  ,1:end-1),n^2,1);
 xe3      = reshape(Xe(1:end-1,2:end  ),n^2,1);
 ye3      = reshape(Ye(1:end-1,2:end  ),n^2,1);
 xe4      = reshape(Xe(2:end  ,2:end  ),n^2,1);
 ye4      = reshape(Ye(2:end  ,2:end  ),n^2,1);


 xl       =  ceil(min(min(xe1/dx,xe3/dx),min(xe2/dx,xe4/dx))+1e-15);
 xr       =  ceil(max(max(xe2/dx,xe4/dx),max(xe1/dx,xe3/dx))-1e-15);
 yd       =  ceil(min(min(ye3/dy,ye4/dy),min(ye1/dy,ye2/dy))+1e-15);
 yu       =  ceil(max(max(ye1/dy,ye2/dy),max(ye3/dy,ye4/dy))-1e-15);



nzind     = [];

for i = 1:n^2

       xlist     = mod(xl(i):xr(i),n);
       ylist     = mod(yd(i):yu(i),n);
       xlist(find(xlist==0))=n;
       ylist(find(ylist==0))=n;


       [Xi,Yi]   = meshgrid(xlist,ylist);
       xi        = reshape(Xi,1,prod(size(Xi)));
       yi        = reshape(Yi,1,prod(size(Yi)));

       nzind     = [ nzind ,((i-1)*n^2)+(yi-1)*n+xi];

end

 nzind   = unique(sort(nzind));
 nnzA    = length(nzind);
 [jind,iind] = ind2sub([n^2,n^2],nzind);

 R1      = sparse(iind,1:nnzA,ones(1,nnzA),n^2,nnzA);
 R2      = sparse(jind,1:nnzA,ones(1,nnzA),n^2,nnzA);

 xpiR    = xpi'*R2;
 R1      = sparse(iind,1:nnzA,xpiR,n^2,nnzA);

 R       = [R1 ; R2];
 b       = [ xpi; ones(n^2,1)    ];
 f       = 0*ones(nnzA,1);  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Av = zeros(nnzA,1);
 for k = 1:nnzA
   i = iind(k);
   j = jind(k);
   px = min(abs(xc(j)-xce(i)),1);
   py = min(abs(yc(j)-yce(i)),1);   

   Av(k) = min(abs( dx-min(px,1-px ) )*abs( dy-min(py,1-py))/dx/dy,1);


 end
  



 Ap = sparse(n^2,n^2);
 Ap(nzind) = Av;
  


 % Cx<d
  C  = [-speye(nnzA), -sparse(ones(nnzA,1)) ;
         speye(nnzA), -sparse(ones(nnzA,1))  ];
  d  =  [ -Av; Av ];
   
  R  = [R sparse(zeros(2*n^2,1))];
  f  = [f ;1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 

 tic
    [xopt,flag] = linprog(f,C,d,R,b,0*ones(nnzA+0,1) ,[]);
 toc
 flag


  xopt           = xopt(1:end-1);

  Avec           = sparse(n^4,1);
  Avec(nzind)    = xopt;
  A              = reshape(Avec,n^2,n^2);
  
 
   plot(Av)
  %Av
  hold on
  plot(A(nzind),'r')
 


  

