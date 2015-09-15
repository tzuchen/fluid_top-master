clear all
%
% fd = 1: solve the problem
% fd = 0: save the problem data
%
% plotstructure = 1: plot it
% plotstructure = 0: don't plot
%
fd = 0;
plotstructure = 1; 
nproc = 4;
den   = 5;
X = createcoord(0,8,32*den,1);
Y = createcoord(0,1,4*den,0);
Z = createcoord(0,1,4*den,0);

cdinfo = createcoords(X,Y,Z);
cdinfo.showtime  = 1;

Lr     = generateL(cdinfo);
Dr     = generateD(cdinfo);
Gr     = Dr';
Pr     = abs(Dr')/2/max(max(Dr));  % Be careful when large size!



% Now comes the BC setting

  BCs{1,1}  = [];               % 0yzu  
  BCs{1,2}  = [];               % 1yzu
  BCs{1,3}  = zeros(10,10);     % x0zu
  BCs{1,4}  = zeros(10,10);     % x1zu
  BCs{1,5}  = zeros(10,10);     % xy0u
  BCs{1,6}  = zeros(10,10);     % xy1u

  BCs{2,1}  = [];               % 0yzv  
  BCs{2,2}  = [];               % 1yzv
  BCs{2,3}  = zeros(10,10);     % x0zv
  BCs{2,4}  = zeros(10,10);     % x1zv
  BCs{2,5}  = zeros(10,10);     % xy0v
  BCs{2,6}  = zeros(10,10);     % xy1v

  BCs{3,1}  = [];               % 0yzw  
  BCs{3,2}  = [];               % 1yzw
  BCs{3,3}  = zeros(10,10);     % x0zw
  BCs{3,4}  = zeros(10,10);     % x1zw
  BCs{3,5}  = zeros(10,10);     % xy0w
  BCs{3,6}  = zeros(10,10);     % xy1w


[lr,dr] = generatelrdr(BCs,cdinfo);


bodyf{1}    =  1;
bodyf{2}    =  0;
bodyf{3}    =  0;

mu = 1;

fr      = generatefr(bodyf,cdinfo);
[A,b]   = generateAb(mu,Lr,Dr,lr,dr,fr,cdinfo);

[p,r]   = domaindecomp([2 2 2],cdinfo);
Ap      = A(p,p);
bp      = b(p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P      = [ Pr ; sparse(cdinfo.nofgrid{cdinfo.dim+1},size(Pr,2))];
Pp     = P(p,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set a preconditioner matrix Ac
% Note: now it is just an identity. When it is not just 
% an identity, DON'T forget the permutation

Ac      = speye(length(p)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set alpha
  clear alphapos
%                xl   xu   yl   yu   zl   zu  
  alphapos{1} = [0.5   1    0.2    0.5    0.2   0.5];
  alphapos{2} = [0.5   1    0.5    0.8    0.5   0.8]; 
  alphapos{3} = [1    1.5   0.2    0.5    0.5   0.8];  
  alphapos{4} = [1    1.5   0.5    0.8    0.2   0.5];
  alphaind = alphaset(alphapos,cdinfo);

 %  alphaijk  = ind2ijk(1:cdinfo.nofgrid{4},4,cdinfo);
 %  alphaind{1} = find(and(abs(alphaijk(1,:)-10-alphaijk(2,:))<2,alphaijk(3,:)<=4));
   

  al2pmat  = 0*al2p(alphaind,cdinfo); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set beta 
% 
%              dir
 betadir{1} = [ 3 ];
 betadir{2} = [ 3 ]; 
 betadir{3} = [ 3 ];  
 betadir{4} = [ 3 ];

 be2almat = be2al(betadir,alphaind,cdinfo);


 be2pmat  = al2pmat*be2almat;
 be2vmat  = Pr*be2pmat;
 be2Amat  = cat(1,be2vmat,sparse(size(A,1)-size(be2vmat,1),size(be2vmat,2)));
 be2Amatp = be2Amat(p,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set objective function. Here is a LINEAR one
%              xl    xu    yl    yu   zl    zu    dir   val
  cpos{1}  =  [0     8    0.3   0.7  0     0.3     2      1  ]; 
  cpos{2}  =  [0     8    0.3   0.7  0.7     1     2     -1  ];
  cpos{3}  =  [8    16    0.3   0.7  0     0.3     2     -1  ]; 
  cpos{4}  =  [8    16    0.3   0.7  0.7     1     2      1  ]; 

  ct = linobjset(cpos,cdinfo);
  c  = 1*cell2mat(ct'); 
  cp = c(p);   % Make it zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set a beta vector 

betalength = size(be2almat,2);
beta       = zeros(betalength,1);
beta(:)    = 20000;
alinp      = find(be2pmat*beta);

%Apb        = diag(sparse(be2Amatp*beta))+Ap;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphavec = be2Amatp*beta;
beinp    = be2pmat*beta;  %design parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the structure

    dy = 0.15;
    [Y0,Z0] = meshgrid(0.15:dy:0.9,0.15:dy:0.9);
    nofline = prod(size(Y0));
    x0 = (cdinfo.coord{1}.l)*ones(1,nofline);
    y0 = reshape(Y0,1,nofline); 
    z0 = reshape(Z0,1,nofline);
    y  = y0;
    z  = z0;
if plotstructure==1
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
   % isosurface(PX,PY,PZ,alphamat);

    axis equal
    axis([0 16 0 1 0 1])


    box on
    grid on
    hold on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now the petsc part

if fd == 1
    disp('Now launching Petsc to solve the problem')
    PetscInitialize(nproc,'fluidsolver');
    socketp = openport();
 
    send(socketp,Ap);
    send(socketp,Ac);
    send(socketp,Pp);
    send(socketp,bp);
    send(socketp,cp);

    %loop here
    for ncount = 1:2
          send(socketp,ncount); % Send the task number
 	  send(socketp,diag(sparse(alphavec)));
   	  fobj  = receive(socketp); % Receive fobj
   	  fobj  = fobj(1);   
   	  xp    = receive(socketp); % Receive xp
   	  x     = xp(r); 
   	  gradf = receive(socketp); % Receive gradf
    end

    PetscFinalize(socketp);


     sol = x;

    % we don't need pressure
    sol(end-cdinfo.nofgrid{4}+1:end)=0;




    [y,z,S,sm] = streamlinemap(y0,z0,cdinfo,sol);
 
    for i=1:size(sm.x,2)
       h = plot3(sm.x(:,i),sm.y(:,i),sm.z(:,i));
       set(h,'color','b','linewidth',1.5)
    end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else %fd = 0;
% write to files
   disp('Now saving the problem data to files')
   R   = sparse(r,1:length(r),ones(length(r),1),length(r),length(r));
   send2file('Apfiled',Ap);
   send2file('Acfiled',Ac);
   send2file('Ppfiled',Pp);
   send2file('bpfiled',bp);
   send2file('cpfiled',cp);
   send2file('alphafiled',diag(sparse(alphavec)));
   send2file('beinpfiled',beinp);
   send2file('Rfiled',R);
   cdinfosend2file('cdinfofiled',cdinfo);

   send2file('x0filed',x0);
   send2file('y0filed',y0);
   send2file('z0filed',z0);

   eval('!chmod +r *')
   eval('!chmod -t *')
   eval('!chmod +w *')

end 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







