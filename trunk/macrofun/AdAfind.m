function [Am,dA]=AdAfind(Y0,Z0,Ye,Ze,Ycor,Zcor,ny,nz,sol,cdinfo,D)

ydim    = cdinfo.coord{2}.r - cdinfo.coord{2}.l;
zdim    = cdinfo.coord{3}.r - cdinfo.coord{3}.l;

dy      = ydim / ny;
dz      = zdim / nz;


y0 = reshape(Y0,ny*nz,1);
z0 = reshape(Z0,ny*nz,1);
y  = reshape(Ye,ny*nz,1);
z  = reshape(Ze,ny*nz,1);

 ye1      = reshape(Ycor(1:end-1,1:end-1),ny*nz,1);
 ze1      = reshape(Zcor(1:end-1,1:end-1),ny*nz,1);
 ye2      = reshape(Ycor(2:end  ,1:end-1),ny*nz,1);
 ze2      = reshape(Zcor(2:end  ,1:end-1),ny*nz,1);
 ye3      = reshape(Ycor(1:end-1,2:end  ),ny*nz,1);
 ze3      = reshape(Zcor(1:end-1,2:end  ),ny*nz,1);
 ye4      = reshape(Ycor(2:end  ,2:end  ),ny*nz,1);
 ze4      = reshape(Zcor(2:end  ,2:end  ),ny*nz,1);

 yl       =  max(min(min(ceil(ye1/dy),ceil(ye3/dy)),min(ceil(ye2/dy),ceil(ye4/dy))),1);
 yr       =  min(max(max(ceil(ye2/dy),ceil(ye4/dy)),max(ceil(ye1/dy),ceil(ye3/dy))),ny);
 zd       =  max(min(min(ceil(ze3/dz),ceil(ze4/dz)),min(ceil(ze1/dz),ceil(ze2/dz))),1);
 zu       =  min(max(max(ceil(ze1/dz),ceil(ze2/dz)),max(ceil(ze3/dz),ceil(ze4/dz))),nz);


nzind     = [];
nofline = ny*nz;


for i = 1:nofline
     if yl(i)==yr(i)
       yl(i)=max(1,yl(i)-1);
       yr(i)=min(ny,yr(i)+1);
     end
     if zd(i)==zu(i)
       zd(i)=max(1,zd(i)-1);
       zu(i)=min(nz,zu(i)+1);
     end
 
     [Yi,Zi]=meshgrid(yl(i):yr(i),zd(i):zu(i));
     yi = reshape(Yi,1,prod(size(Yi)));
     zi = reshape(Zi,1,prod(size(Yi)));
     nzind         = [ nzind ,((i-1)*ny*nz)+(yi-1)*nz+zi];

end
nzind = unique(sort(nzind));
[iind,jind] = ind2sub([ny*nz,ny*nz],nzind);
nnzA = length(nzind);


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
 VX0      = 1000*interp2(cY,cZ,vx0a,Y0,Z0);


 b       = [ reshape(VX0,ny*nz,1); reshape(VX0,ny*nz,1) ];
 R1      = sparse(iind,1:nnzA,ones(1,nnzA),ny*nz,nnzA);
 R2      = sparse(jind,1:nnzA,ones(1,nnzA),ny*nz,nnzA);
 R       =  [ R1 ;  R2];  
 f       = 0*ones(nnzA,1);


 Mv = zeros(nnzA,1);
 for k = 1:nnzA
   i = iind(k);
   j = jind(k);
   Mv(k) = abs( abs(y0(i)-y(j)) - dy )*abs( abs(z0(i)-z(j)) - dz );
   Mv(k) = Mv(k)*VX0(i)*ny*nz;
 end
   
 Mp = sparse(ny*nz,ny*nz);
 Mp(nzind) = Mv;

 % Cx<d
  C  = [-speye(nnzA), -sparse(ones(nnzA,1)) ;
         speye(nnzA), -sparse(ones(nnzA,1))  ];
  d  = [ -Mv; Mv ];
   

  R  = [R sparse(zeros(2*ny*nz,1))];
  f  = [f ;1];

tic
  [f,fv,flag] =linprog(f,C,d,R,b,0*ones(nnzA,1),[]);
toc
  flag



  Avec = sparse(ny^2*nz^2,1);
  Avec(nzind) = f(1:nnzA);
  Am          = reshape(Avec,ny*nz,ny*nz);

  for i = 1:ny*nz
    
    Am(i,:) = Am(i,:)/ sum(Am(i,:));
  end

finddA = 0;
dA = [];
if finddA == 1;
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

  %yave    =(yivec+yjvec)/2;
  %zave    =(zivec+zjvec)/2;
  %vave    = interp2(cY,cZ,vx0a,yave,zave);

  wdydv   = (dz - abs(zdiff)).*sign(ydiff).*dlambda';
  wdzdv   = (dy - abs(ydiff)).*sign(zdiff).*dlambda';

  S = sparse(indj,1:nnzA,1);
  dA = (S*wdydv)' *D.dydv + (S*wdzdv)' *D.dzdv; 
end

