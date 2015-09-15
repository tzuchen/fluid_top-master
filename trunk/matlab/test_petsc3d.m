mypathdef

% main routine
 close all
 loaddata = 0;
if loaddata==0
	n              = 30;
 	matratio       = 0.6;
 	minfo.period   = 1;
 	minfo.alphacon = 0;
 	minfo.iternum  = 10;
 	minfo.BC       = @BC_parabolic3d;
 	minfo.solver   = 'symmlq';
	%  
 	meshdata = generatemesh3d(n,minfo);
 	meshdata = BCsetting(meshdata,minfo);
 	meshdata = ObjectiveSetting(meshdata,minfo);
 	meshdata = generateoperator3df(meshdata,minfo);

 	meshdata = petscAbgenerate_3d(n,meshdata,minfo);
        save meshdata
else
        load meshdata40

end
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.fnamealpha = 'alphavec.mat';
opt.rtol       =  1e-5;
opt.ngrid      =  n;
opt.ksp_type    =  'cg'
err = PetscInitialize(4,opt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here let's generate a shaped tube
   betamat1 = eye(n)+diag(ones(n-1,1,1),1)+diag(ones(n-1,1),-1);  
   betamat2 = betamat1(n:-1:1,:);
   betavec  = zeros(2,2*n^2);
   betavec(1,:) = [reshape(betamat1,n^2,1);reshape(betamat2,n^2,1)]';
   betavec(2,:) = [reshape(betamat2,n^2,1);reshape(betamat1,n^2,1)]';


   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
meshdata.fnamealpha = 'alphavec.mat';meshdata.port  = 0;
% solve this 2 times
for count = 1:2
    meshdata.betavec  = 8000*betavec(mod(count,2)+1,:);
    meshdata = beta2alpha2(meshdata);
    meshdata = diffsolverpetsc3d(meshdata,1,minfo);
end
meshdata.alphamat= reshape(meshdata.alphavec,n,n,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
port = meshdata.port;
err = PetscFinalize(port);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plotffield3d(meshdata,0);


