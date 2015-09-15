% diffsolver
function  meshstruc = diffsolver(meshstruc,alphagain,info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lr       = meshstruc.Lr;
lr       = meshstruc.lr;
Gr       = meshstruc.Gr;
Dr       = meshstruc.Dr;
Pr       = meshstruc.Pr;
dr       = meshstruc.dr;
fr       = meshstruc.fr;
alphavec = meshstruc.alphavec;
mu       = meshstruc.mu;
n        = meshstruc.n;
nr       = meshstruc.nr;
nur      = meshstruc.nur;
nvr      = meshstruc.nvr;
nwr      = meshstruc.nwr;
np       = meshstruc.np;
ngd      = meshstruc.ngd;
npt      = meshstruc.npt;
ugd      = meshstruc.ugd;
vgd      = meshstruc.vgd;
wgd      = meshstruc.wgd;
pgd      = meshstruc.pgd;
c        = meshstruc.c;
uindexmat  = meshstruc.uindexmat;
vindexmat  = meshstruc.vindexmat;
windexmat  = meshstruc.windexmat;
pindexmat  = meshstruc.pindexmat;
usollist = meshstruc.usollist;
vsollist = meshstruc.vsollist;
wsollist = meshstruc.wsollist;
sollist  = meshstruc.sollist;
sol      = meshstruc.sol;
adjvec   = meshstruc.adjvec;
objtype  = meshstruc.objtype;
%solver   = info.solver;
port     = meshstruc.port;
r        = meshstruc.r;
p        = meshstruc.p;
fnamealpha = meshstruc.fnamealpha;
alphalist = meshstruc.alphalist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some useful indices
% nr = nur+nvr+nwr
global ind
ind.ur    = 1         : nur;
ind.vr    = nur+1     : nur+nvr;
ind.wr    = nur+nvr+1 : nur+nvr+nwr;

ind.velr     = 1         :nr;
ind.pr       = nr+1      :nr+np;
%ind.alphar   = nr+np+1   :nr+2*np;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 cr        =  [ c(sollist); zeros(np,1) ];
 crp       =  cr(p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if port ==0
   port = openport(5050);    
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

send(port, alphavec);
%send(crp);
%send(alphalist);

%PetscBinaryWrite(fnamealpha,alphavec);
%PetscBinaryWrite('cr.mat',crp);
%PetscBinaryWrite('alphalist.mat',alphalist);

%taskflag = 1;
%save('taskflag.mat','taskflag');  
 

   
    fobj  = receive(port);
    x     = receive(port);
    gradf = receive(port);
    
    x = x(r);
    dfdalpha = sparse(np,1);
    dfdalpha(alphalist) = -gradf;
    disp('One set of f and Df is solved!')
    
        sol    = [x ;0];
        obj    = fobj(1);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   ur       = sol(ind.ur);
   vr       = sol(ind.vr);
   wr       = sol(ind.wr);
   pr       = sol(ind.pr);
  % alphar   = sol(ind.alphar);      
   
   u        = zeros(n, 1);
   v        = zeros(n, 1);   
   w        = zeros(n, 1);
   p        = zeros(np,1);
    
   u(usollist) = ur;
   v(vsollist) = vr;
   w(wsollist) = wr;   
   p           = pr;
   %obj         = caug(1:end-1)'*sol;
   
   for ncount = usollist
       ugd{ncount}.u = u(ncount);           
   end
   for ncount = vsollist
       vgd{ncount}.v = v(ncount);
   end
   for ncount = wsollist
       wgd{ncount}.w = w(ncount);
   end
      
   for npcount = 1:np
       pgd{npcount}.p  =  p(npcount);
   end
   
   if info.period == 1
   % period condition
        for ncount = 1:npt*ngd*ngd
              i = ugd{ncount}.i;
              j = ugd{ncount}.j;
              k = ugd{ncount}.k;
              if i == npt
                 ugd{ncount}.u = ugd{uindexmat(1  ,j , k)}.u;
              end   
        end
   end
   
   umat  = zeros(npt,ngd,ngd);
   vmat  = zeros(ngd,npt,ngd);
   wmat  = zeros(ngd,ngd,npt);
   pmat  = zeros(ngd,ngd,ngd);
   
   for ncount = 1 : n
        umat(ugd{ncount}.i,ugd{ncount}.j,ugd{ncount}.k) = ugd{ncount}.u;
        vmat(vgd{ncount}.i,vgd{ncount}.j,vgd{ncount}.k) = vgd{ncount}.v;
        wmat(wgd{ncount}.i,wgd{ncount}.j,wgd{ncount}.k) = wgd{ncount}.w;        
   end
   
   for npcount = 1:np
      pmat(pgd{npcount}.i,pgd{npcount}.j,pgd{npcount}.k)  = pgd{npcount}.p;   
   end
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for npcount = 1 : np
    i = pgd{npcount}.i;
    j = pgd{npcount}.j;
    k = pgd{npcount}.k;
    pgd{npcount}.u = 1/2*(ugd{uindexmat(i+1,j  ,k  )    }.u + ugd{uindexmat(i,j,k)    }.u);
    pgd{npcount}.v = 1/2*(vgd{vindexmat(i  ,j+1,k  )-n  }.v + vgd{vindexmat(i,j,k)-n  }.v);
    pgd{npcount}.w = 1/2*(wgd{windexmat(i  ,j  ,k+1)-2*n}.w + wgd{windexmat(i,j,k)-2*n}.w);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate Energy
  %energy = -[ur;vr;wr]'*(mu*Lr-1*diag(Pr*alphavec))*[ur;vr;wr];
   energy = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meshstruc.pgd   = pgd;
meshstruc.ugd   = ugd;
meshstruc.vgd   = vgd;
meshstruc.wgd   = wgd;
meshstruc.umat  = umat;
meshstruc.vmat  = vmat;
meshstruc.wmat  = wmat;
meshstruc.pmat  = pmat;
meshstruc.ur   = ur;
meshstruc.vr   = vr;
meshstruc.wr   = wr;
meshstruc. cr   = cr; 
meshstruc.energy= energy;
meshstruc.sol   = sol;
meshstruc.obj   = obj;
meshstruc.port   = port;
meshstruc.dfdalpha  = dfdalpha;
