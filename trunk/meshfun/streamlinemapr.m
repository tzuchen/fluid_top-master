function [y,z,D,sm] = streamlinemapr(y0,z0,cdinfo,sol)
% given a set of points y0,z0 and x0=0, this function calculate the 
% streamlines and return the end points y and z.
%
% Notice:
%  cdinfo and sol are required !
%  This function calls petsc routine "streamliner"
%
% Tzu-Chen Liang 6/12/2006
  
   maxnum  = 5000;    
   w = size(y0);
   y0 = reshape(y0,prod(w),1);
   z0 = reshape(z0,prod(w),1);
   totalnofline = prod(w);
   y = y0;
   z = z0;

   ngroup = ceil((totalnofline-1)/maxnum);
   ind = 1;

   for i = 1:ngroup
 
    if  i<ngroup
     nofline = maxnum;
    else
     nofline = totalnofline-ind+1;
    end

     yi = y0(ind:ind+nofline-1);
     zi = z0(ind:ind+nofline-1);



   PetscInitialize(5,'streamliner');
   socketp = openport();
   
   cdinfosend(socketp,cdinfo);
   send(socketp,sol);



 
   
   xi = (cdinfo.coord{1}.l)*ones(1,nofline);
 
   
   send(socketp,xi);
   send(socketp,yi);
   send(socketp,zi);

   xe      = receive(socketp);
   ye      = receive(socketp);
   ze      = receive(socketp);   
   
   %sx      = receive(socketp);
   %sy      = receive(socketp);
   %sz      = receive(socketp);

   D.dydv = SparseReceive(socketp);
   D.dydv = D.dydv(:,1:end-24);
   D.dzdv = SparseReceive(socketp);
   D.dzdv = D.dzdv(:,1:end-24);

   n = 500;
   %nofline = length(y0);


   %sm.x = reshape(sx,n,nofline);
   %sm.y = reshape(sy,n,nofline);
   %sm.z = reshape(sz,n,nofline);
   sm = [];

   y(ind:ind+nofline-1)= ye;
   z(ind:ind+nofline-1)= ze;
   ind = ind+nofline;
 

   %PetscFinalize(socketp);
   closeport(socketp);

 end


   y = reshape(y,w);
   z = reshape(z,w);
