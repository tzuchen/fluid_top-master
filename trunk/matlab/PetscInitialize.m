function err = PetscInitialize(np,option)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Set options
% default value 
  rtol     = 1e-10;
  ksp_type = 'symmlq';
  ngrid    = 20;
  fnameA    = 'A20.mat';
  fnameb    = 'b20.mat'; 
  fnameAp   = 'Ap20.mat';
  fnamePr   = 'Pr20.mat';
  fnamealpha = 'alphavec20.mat';
  fparam     = 'param20';
  
 if nargin == 2
    if isfield(option,'rtol')
        rtol   = option.rtol;
    end
    if isfield(option,'ksp_type')
        ksp_type   = option.ksp_type;
    end
    if isfield(option,'ngrid')
        ngrid      = option.ngrid;
        fnameA      = ['A',num2str(ngrid),'.mat'];
        fnameAp     = ['Ap',num2str(ngrid),'.mat'];
        fnameb      = ['b',num2str(ngrid),'.mat'];  
        fnamePr     = ['Pr',num2str(ngrid),'.mat'];
        fnamealpha     = ['alphavec',num2str(ngrid),'.mat'];
        fparam         = ['param',num2str(ngrid),'.mat'];
    end    
    
    if isfield(option,'fnamealpha')
     fnamealpha     = option.fnamealpha;
    end
    
 end
    
    load(fparam)
    
 fparam
    opt{1}  = [' -rtol ', num2str(rtol) ]; 
    opt{2}  = [' -ksp_type ', ksp_type ]; 
    opt{3}  = [' -ngrid ', num2str(ngrid) ];
    opt{4}  = [' -fnameA ', fnameA ];
    opt{5}  = [' -fnameAp ', fnameAp ];
    opt{6}  = [' -fnameb ', fnameb ];
    opt{7}  = [' -fnamePr ', fnamePr ];
    opt{8}  = [' -fnamealpha ', fnamealpha ];
    opt{9}  = [' -ngrid ' , num2str(param(1))];
    opt{10} = [' -neq ' , num2str(param(2))];   
    opt{11} = [' -na '  , num2str(param(3))];
    opt{12} = [' -n3 ' ,  num2str(param(4))];
    opt{13} = [' -nalpha ' ,  num2str(param(5))];    
    opt{14} = [' -log_summary '];
    optset = [opt{:}];
opt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist('taskflag.mat','file')
   delete taskflag.mat
end
if exist('endflag.mat','file')
   delete endflag.mat
end

err = launch(['./mysolver' ' -viewer_socket_machine ' getenv('HOST') optset],np);
if err==0
   disp(sprintf('Petsc solver launched successfully!'))
else
   disp(sprintf('Failed to launch Petsc solver! Please stop the program.'))
end



