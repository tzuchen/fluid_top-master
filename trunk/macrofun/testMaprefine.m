%  Test routine for Maprefine.c
%  
%
% 
%  Tzu-Chen Liang 9-28-2006

nproc = 3;
n     = 500;
iter  = 5;

option.withMatlab = 1;
PetscInitialize(nproc,'Maprefine2',option);
socketp = openport();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEND part %
 send(socketp,n);
 send(socketp,iter);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECEIVE part %

x0      = receive(socketp);
%y0      = receive(socketp);
mixnorm =  receive(socketp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



mixnorm = mixnorm(1:iter);

imagesc(reshape(x0,n/2,n))
mixnorm


PetscFinalize(socketp);

