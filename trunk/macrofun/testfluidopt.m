%  Test routine for fluidopt.c
%  
%
% 
%  Tzu-Chen Liang 9-28-2006

nproc = 4;


option.withMatlab = 1;
option.rtol       = 1e-10;
option.abstol     = 1e-10;
PetscInitialize(nproc,'fluidopt',option);
socketp = openport();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEND part %
% Nothing to send

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECEIVE part %

for i = 1:1
disp(i)
sol  = receive(socketp);
A = SparseReceive(socketp,2500,2500);
y  = receive(socketp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




PetscFinalize(socketp);
