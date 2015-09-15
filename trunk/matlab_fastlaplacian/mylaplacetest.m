close all
n             = 4;
info.matratio = 0.6;
info.period   = 1;
info.alphacon = 0;
info.alphamin = 0;
info.alphamax = 8000;

info.BC       = @BC_parabolic3d;
info.solver   = 'symmlq';


tic
meshdata = generatemesh3d(n,info);
toc
tic
meshdata = BCsetting(meshdata,info);
toc
tic
meshdata = generateoperator3df(meshdata,info);
toc

tic
meshdata = petscAbgenerate_3d(n,meshdata,info);
toc
