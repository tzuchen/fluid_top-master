function petsctest3rfun2(n,filename,betainitial)
%%%%%%%%%%%%%%%%%%%%%
% path setting (this is not good!)
mypathdef
%%%%%%%%%%%%%%%%%%%%%

global meshdata
global info

disp(sprintf('Initializing Petsc/Snopt version optimizer.'))


    %n             = 30;
    info.matratio = 0.6;
    info.period   = 1;
    info.alphacon = 0;
    info.alphamin = 0;
    info.alphamax = 8000;

    info.BC       = @BC_parabolic3d;
    info.solver   = 'symmlq';
if nargin == 1 % no precalculated data 
    meshdata = generatemesh3d(n,info);
    meshdata = BCsetting(meshdata,info);
    meshdata = ObjectiveSetting(meshdata,info);
    meshdata = generateoperator3ds(meshdata,info);
    meshdata = petscAbgenerate_3d(n,meshdata,info);
    disp(sprintf('Meshdata is generated!'))
    
    filename =  ['meshdata',num2str(n)];
    save(filename,'meshdata');
    disp(sprintf('Meshdata is saved (for next time use)!'))
    
else
    load(filename)     
    disp(sprintf('Meshdata is loaded!'))
end
    nalpha    = meshdata.nalpha;
    pgd       = meshdata.pgd;
    alphalist = meshdata.alphalist;
    alphavec  = meshdata.alphavec;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial setting
x      = zeros(nalpha,1);
if info.alphacon ==1
    meshdata = alphasetting(meshdata,info);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.fnamealpha = 'myalpha.mat';
opt.rtol       =  1e-10;
opt.ngrid      =  n;
err = PetscInitialize(4,opt);
meshdata.fnamealpha = 'myalpha.mat';
meshdata.port  = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  set initial beta
if nargin > 2
    disp(sprintf('Load initial betavec!'))
    load(betainitial)
    meshdata.betavec = backup.betavec; 
    meshdata.nbeta   = length(meshdata.betavec);
    meshdata         = dalphadbeta2(meshdata);
else
    disp(sprintf('Generate a zero initial betavec!'))
    meshdata.nbeta   = 2*n^2;
    meshdata.betavec = 0*ones(meshdata.nbeta,1);
    meshdata         = dalphadbeta2(meshdata);
end
x        = meshdata.betavec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meshdata  =  ObjectiveSetting(meshdata,info);
meshdata = diffsolverpetsc3d(meshdata,1,info);
meshdata.port
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % optimization part
[x,F,inform] = snoptflow3r2(meshdata,x);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
port = meshdata.port
err = PetscFinalize(port);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alphavec  = meshdata.alphavec;

alphamat  = zeros(n,n,n);

 for pcount = 1:meshdata.np
     pgd{pcount}.alpha                        = alphavec(pcount);
     alphamat(pgd{pcount}.i,pgd{pcount}.j,pgd{pcount}.k)    = alphavec(pcount);
 end
% 
%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meshdata.alphavec = alphavec;
meshdata.alphamat = alphamat;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%meshdata = diffsolver3d(meshdata,1,info);
%meshdata = diffsolverpetsc3d(meshdata,1,info);
%plotffield3d(meshdata,0);

save meshdatasnopt

