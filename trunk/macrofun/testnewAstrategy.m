function [M,Am,Mp,dA,mdot,tempind] = testnewAstrategy(n,sol,cdinfo)

ydim    = cdinfo.coord{2}.r - cdinfo.coord{2}.l;
zdim    = cdinfo.coord{3}.r - cdinfo.coord{3}.l;

dy      = ydim / n;
dz      = zdim / n;

ypt     = dy/2:dy:ydim-dy/2;
zpt     = dz/2:dz:zdim-dz/2;

[Y0,Z0] = meshgrid(ypt,zpt);
nofline = prod(size(Y0)); 

x0      = (cdinfo.coord{1}.l)*ones(1,nofline);
y0      = reshape(Y0,1,nofline); 
z0      = reshape(Z0,1,nofline);


[y,z,D] = streamlinemap(y0,z0,cdinfo,sol);

ye       = 0:dy:1;
ze       = 0:dz:1;
[Ye,Ze]  = meshgrid(ye,ze);
ye       = reshape(Ye,(n+1)^2,1);
ze       = reshape(Ze,(n+1)^2,1);

[yea,zea,spm] = streamlinemap(ye,ze,cdinfo,sol);
Ye       = reshape(yea,n+1,n+1);
Ze       = reshape(zea,n+1,n+1);




 ye1      = reshape(Ye(1:end-1,1:end-1),n^2,1);
 ze1      = reshape(Ze(1:end-1,1:end-1),n^2,1);
 ye2      = reshape(Ye(2:end  ,1:end-1),n^2,1);
 ze2      = reshape(Ze(2:end  ,1:end-1),n^2,1);
 ye3      = reshape(Ye(1:end-1,2:end  ),n^2,1);
 ze3      = reshape(Ze(1:end-1,2:end  ),n^2,1);
 ye4      = reshape(Ye(2:end  ,2:end  ),n^2,1);
 ze4      = reshape(Ze(2:end  ,2:end  ),n^2,1);

 yl       =  max(min(min(ceil(ye1/dy),ceil(ye3/dy)),min(ceil(ye2/dy),ceil(ye4/dy))),1);
 yr       =  min(max(max(ceil(ye2/dy),ceil(ye4/dy)),max(ceil(ye1/dy),ceil(ye3/dy))),n);
 zd       =  max(min(min(ceil(ze3/dz),ceil(ze4/dz)),min(ceil(ze1/dz),ceil(ze2/dz))),1);
 zu       =  min(max(max(ceil(ze1/dz),ceil(ze2/dz)),max(ceil(ze3/dz),ceil(ze4/dz))),n);

nzind     = [];
meanind   = [];
for i = 1:nofline
     [Yi,Zi]=meshgrid(yl(i):yr(i),zd(i):zu(i));
     yi = reshape(Yi,1,prod(size(Yi)));
     zi = reshape(Zi,1,prod(size(Yi)));
     nzind         = [ nzind ,((i-1)*n^2)+(yi-1)*n+zi];
 
end


nzind       = unique(sort(nzind));
nnzA        = length(nzind);
[iind,jind] = ind2sub([n^2,n^2],nzind);

     ymean   = ceil(y/dy);
     zmean   = ceil(z/dz);
     tempind = sub2ind([n,n],ymean,zmean);
for i = 1:nofline
 meanind         = [ meanind ,((i-1)*n^2)+(ymean(i)-1)*n+zmean(i)];
end


[nouse,meanind]=intersect(nzind,meanind);


 vxspam   = cdinfo.gridspam{1};
 vx       = reshape(sol(1:cdinfo.nofgrid{1}),vxspam(1),vxspam(2),vxspam(3));
 vx0 	  = squeeze(vx(1,:,:)+vx(end,:,:))/2;
 vx0a     = zeros(vxspam(2)+2,vxspam(3)+2);
 vx0a(2:end-1,2:end-1) = vx0;
 cdy 	  = ydim / vxspam(2);
 cdz 	  = zdim / vxspam(3);
 cy 	  = -cdy/2:cdy:ydim+cdy/2;
 cz       = -cdz/2:cdz:zdim+cdz/2;
 [cY,cZ]  = meshgrid(cy,cz);
 VX0      = 1*interp2(cY,cZ,vx0a,Y0,Z0);
 VX0 = VX0./sum(sum(VX0));
 
 
 b       = [ reshape(VX0,n^2,1); reshape(VX0,n^2,1) ];
 R1      = sparse(iind,1:nnzA,ones(1,nnzA),n^2,nnzA);
 R2      = sparse(jind,1:nnzA,ones(1,nnzA),n^2,nnzA);
 R       =  [ R1 ;  R2];  
 f       = 0*ones(nnzA,1);
 
 Mv = zeros(nnzA,1);
 for k = 1:nnzA
   i = iind(k);
   j = jind(k);
   Mv(k) = abs( abs(y0(i)-y(j)) - dy )*abs( abs(z0(i)-z(j)) - dz );
   Mv(k) = Mv(k)*VX0(i)*n^2;
 end
   
 Mp = sparse(n^2,n^2);
 Mp(nzind) = Mv;

 % Cx<d
  C  = [-speye(nnzA), -sparse(ones(nnzA,1)) ;
         speye(nnzA), -sparse(ones(nnzA,1))  ];
  d  =  [ -Mv; Mv ];
   

  R  = [R sparse(zeros(2*n^2,1))];
  f  = [f ;1];



    
tic
      [xopt,fv,flag] =linprog(f,C,d,R,b,zeros(nnzA,1),[]);
 
  %   [xopt,fv,flag] =linprog(f,[],[],R,b,zeros(nnzA,1),[]);
 
toc
  flag

  
  mdot  = R1*xopt(1:end-1);

  xopt = xopt(1:end-1);


  Avec = sparse(n^4,1);
  Avec(nzind) = xopt(1:nnzA);
  Am          = reshape(Avec,n^2,n^2);
  M           = Am;

  for i = 1:n^2
    Am(i,:) = Am(i,:)/ sum(Am(i,:));
  end

finddA=0;
dA = [];
if finddA ==1
  OPTS.disp = 0;
  [Va,Da]       = eigs(Am',6,'LM',OPTS);
  [Dr,I]        = sort(max(abs(Da)));
  Im            = I(end-2);
  mu            = Da(I(end-2),I(end-2));
  dlambda       = real(Va(:,Im))*real(Va(:,Im))';
  dlambda       = dlambda(nzind);


[indi,indj] = find(Am);

yivec   = y(indi);
zivec   = z(indi);
yjvec   = y0(indj)';
zjvec   = z0(indj)';


zdiff   = zivec-zjvec;
ydiff   = yivec-yjvec;

yave    =(yivec+yjvec)/2;
zave    =(zivec+zjvec)/2;
vave    = interp2(cY,cZ,vx0a,yave,zave);



wdydv   = (dz - abs(zdiff)).*sign(ydiff).*dlambda'.*vave;
wdzdv   = (dy - abs(ydiff)).*sign(zdiff).*dlambda'.*vave;



S = sparse(indj,1:nnzA,1);
dA = (S*wdydv)' *D.dydv + (S*wdzdv)' *D.dzdv; 

end

















