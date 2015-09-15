function meshstruc = petscAbgenerate_3d(n,meshstruc,kinfo)

np  = meshstruc.np;
meshstruc.alphavec = 0*ones(np,1); 
meshstruc = diffeq3d(meshstruc,1,kinfo);
 
      
A  = meshstruc.A;
b  = meshstruc.b;
cr = meshstruc.cr(1:end-1);
Pr = meshstruc.Pr;

neq = length(b);
na  = size(Pr,1);
n3  = size(Pr,2);
nalpha = length(meshstruc.alphalist);

Pz  = sparse(neq-na,n3);

PrL = [Pr;Pz];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Re-order A
%p = symrcm(A);

p = permutegenerate(n);
A = A(p,p);
b = b(p);
cr = cr(p);
PrL = PrL(p,:);

meshstruc.p = p;
r = zeros(1,neq);
r(p) = 1:neq;
meshstruc.r = r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ap = sparse(diag(sparse(ones(neq,1))));
param = [meshstruc.ngd ;
         neq ;
         na ;
         meshstruc.ngd^3;
         nalpha];





PetscBinaryWrite(['b',num2str(n),'.mat'],b);
PetscBinaryWrite(['A',num2str(n),'.mat'],A);
PetscBinaryWrite(['Ap',num2str(n),'.mat'],Ap);
PetscBinaryWrite(['Pr',num2str(n),'.mat'],PrL);
%PetscBinaryWrite(['cr',num2str(n),'.mat'],cr);
PetscBinaryWrite(['cr.mat'],cr);
PetscBinaryWrite('alphalist.mat',meshstruc.alphalist);

pfname = ['param',num2str(n)];
save(pfname,'param');

fname = ['meshstruc',num2str(n)];
%save(fname,'meshstruc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

