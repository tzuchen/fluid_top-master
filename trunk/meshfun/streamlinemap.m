function [y,z,D,sm] = streamlinemap(y0,z0,cdinfo,sol)
% given a set of points y0,z0 and x0=0, this function calculate the 
% streamlines and return the end points y and z.
%
% Notice:
%  cdinfo and sol are required !
%  This function calls petsc routine "streamline"
%
% Tzu-Chen Liang 1/12/2006

   


   PetscInitialize(1,'streamline');
   socketp = openport();
   
   cdinfosend(socketp,cdinfo);
   send(socketp,sol);

 
   nofline = length(y0); 
   x0 = (cdinfo.coord{1}.l)*ones(1,nofline);
 
   
   send(socketp,x0);
   send(socketp,y0);
   send(socketp,z0);

   x      = receive(socketp);
   y      = receive(socketp);
   z      = receive(socketp);   
   
   sx      = receive(socketp);
   sy      = receive(socketp);
   sz      = receive(socketp);

   D.dydv = SparseReceive(socketp);
   D.dydv = D.dydv(:,1:end-24);
   D.dzdv = SparseReceive(socketp);
   D.dzdv = D.dzdv(:,1:end-24);

   n = 500;
   nofline = length(y0);


   sm.x = reshape(sx,n,nofline);
   sm.y = reshape(sy,n,nofline);
   sm.z = reshape(sz,n,nofline);



   %PetscFinalize(socketp);
   closeport(socketp);
