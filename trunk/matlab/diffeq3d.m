% diffsolver
function  meshstruc = diffeq3d(meshstruc,alphagain,info)
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

np       = meshstruc.np;
ngd      = meshstruc.ngd;
npt      = meshstruc.npt;
cr       = meshstruc.cr;

adjvec   = meshstruc.adjvec;
objtype  = meshstruc.objtype;

solver   = info.solver;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   A        = -[ mu*Lr-alphagain*diag(sparse(Pr*sparse(alphavec)))            -Gr    ; 
                   -Gr'                                           sparse(np,np)  ];
            
   b        = -[ (-mu*lr - fr) ;
                     dr      ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for later use to find dc/dalpha
    caug = zeros(size(A,1),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch objtype 
        case 'linear'
             caug(1:length(cr)) = cr;
        case 'energy'  
             % do nothing
    end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        A = A(1:end-1,1:end-1);
        b = b(1:end-1);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meshstruc.A     = A;
meshstruc.b     = b;
