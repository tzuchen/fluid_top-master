% Now I want a large problem!

%clear all
savefile = 1;

nproc       =  3;
den         =  5;
bodyf{1}    =  0.1;
bodyf{2}    =  0;
bodyf{3}    =  0;
mu          =  1;
ny          = 50;
nz          = ny;

[Ap,Ac,Pp,Pr,bp,p2Amatp,p,r,cdinfo] = problemgen(nproc,den,bodyf,mu);


   alphaijk  = ind2ijk(1:cdinfo.nofgrid{4},4,cdinfo);
  % alphaind{1} = find(and(abs(alphaijk(1,:)-10-alphaijk(2,:))<2,alphaijk(3,:)<=4));
   alphaind{1} = find(and(mod((alphaijk(1,:)-1.5*alphaijk(2,:)-10),10)<3,alphaijk(3,:)<=4));
   alphaind{2} = find(and(mod((alphaijk(1,:)+1.5*alphaijk(2,:)-10),10)<3,alphaijk(3,:)<=4));
  
       ind1 = find(or(and(alphaijk(1,:)<=200,alphaijk(2,:)<=30),and(alphaijk(1,:)>=200,alphaijk(2,:)<=10)));
       ind2 = find(or(and(alphaijk(1,:)<=200,alphaijk(2,:)>=30),and(alphaijk(1,:)>=200,alphaijk(2,:)>=10)));
   alphaind{1} = intersect(alphaind{1},ind1);
   alphaind{2} = intersect(alphaind{2},ind2);

   al2pmat  = al2p(alphaind,cdinfo); 
   alpha       = 1e6*ones(length(cell2mat(alphaind)),1);
   alphainp    = al2pmat*alpha;
   alphavec    = p2Amatp*alphainp; 
   

cp  = bp*0;
n   = size(Ap,1);
R   = sparse(r,1:n,ones(n,1),n,n);

if savefile==1
   send2file('Apfiled',Ap);
   send2file('Acfiled',Ac);
   send2file('Ppfiled',Pp);
   send2file('bpfiled',bp);
   send2file('cpfiled',cp);
   send2file('alphafiled',diag(sparse(alphavec)));
   send2file('Rfiled',R);

   cdinfosend2file('cdinfofiled',cdinfo);
   eval('!chmod +r *')
   eval('!chmod -t *')
end


%[y0,z0]=meshgrid(0.05:0.05:0.95,0.05:0.05:0.95);
[y0,z0]=meshgrid(0.1:0.1:0.9,0.1:0.1:0.9);
 x0    = zeros(size(z0));
send2file('x0filed',x0);
send2file('y0filed',y0);
send2file('z0filed',z0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 betadir{1} = [ 3 ];
 betadir{2} = [ 3 ];
 be2almat = be2al(betadir,alphaind,cdinfo);
 be2pmat  = al2pmat*be2almat;
 be2vmat  = Pr*be2pmat;
 be2Amat  = cat(1,be2vmat,sparse(size(Ap,1)-size(be2vmat,1),size(be2vmat,2)));
 be2Amatp = be2Amat(p,:);
 betalength = size(be2almat,2);
 beta       = zeros(betalength,1);
 beta(:)    = 20000;
 alinp      = find(be2pmat*beta);

 close all
 pgrid = cdinfo.gridspam{4}
 alphamat        = zeros(pgrid);
 alphamat(alinp) = 1;
 alphamat = shiftdim(alphamat,2);
 alphamat = permute(alphamat,[3 2 1]);


 px  = cdinfo.coord{1}.l+0.5*cdinfo.coord{1}.gsize:cdinfo.coord{1}.gsize:cdinfo.coord{1}.r-0.5*cdinfo.coord{1}.gsize;
 py  = cdinfo.coord{2}.l+0.5*cdinfo.coord{2}.gsize:cdinfo.coord{2}.gsize:cdinfo.coord{2}.r-0.5*cdinfo.coord{2}.gsize;
 pz  = cdinfo.coord{3}.l+0.5*cdinfo.coord{3}.gsize:cdinfo.coord{3}.gsize:cdinfo.coord{3}.r-0.5*cdinfo.coord{3}.gsize;
 [PX,PY,PZ]=meshgrid(px,py,pz);
 isosurface(PX,PY,PZ,alphamat);

 camlight; lighting phong
 hold on

 if exist('output1')==1
   
    A = reshape(output1,500,361,3);
   for i = 2:38:361
    plot3(A(:,i,1),A(:,i,2),A(:,i,3));
    hold on
   end



   for i = 4:38:361
    plot3(A(:,i,1),A(:,i,2),A(:,i,3));
    hold on
   end

 end


 box on
 grid on
 axis([0 10 0 1 0 1])
 axis equal
 axis([0 10 0 1 0 1])


%%%%%%%%%%%%%%%%%%

 for i = 1:1:361
    plot3(A(:,i,1),A(:,i,2),A(:,i,3));
    hold on
   end


