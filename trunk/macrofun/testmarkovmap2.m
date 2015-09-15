function [A] = testmarkovmap2(m,n,cdinfo,sol)
% Test function to find A and dA
%
%
% Tzu-Chen Liang  2-16-2006


   PetscInitialize(4,'testmarkovmap2');
   socketp = openport();
   
   % send Cdinfo and sol
   cdinfosend(socketp,cdinfo);
   send(socketp,sol);
   send(socketp,m);
   send(socketp,n);

   OPTS.disp = 0;

   % just for it to function correctly
   [A]         = SparseReceive(socketp,m*n,m*n); 

 
   closeport(socketp);
