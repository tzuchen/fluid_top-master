function [A,dA,mu,D] = testmarkovmap(m,n,cdinfo,sol)
% Test function to find A and dA
%
%
% Tzu-Chen Liang  2-16-2006


   PetscInitialize(5,'testmarkovmap');
   socketp = openport();
   
   % send Cdinfo and sol
   cdinfosend(socketp,cdinfo);
   send(socketp,sol);
   send(socketp,m);
   send(socketp,n);

   OPTS.disp = 0;

   % just for it to function correctly
   [A]         = SparseReceive(socketp,m*n,m*n); 

   [V,D]       = eigs(A',6,'LM',OPTS);
   [Dr,I]      = sort(max(abs(D)));
   Im          = I(end-1);
   mu          = D(I(end-1),I(end-1));
   nzInd       = find(A);
   dlambda     = real(V(:,Im))*real(V(:,Im))';
   dlambda     = dlambda(nzInd);

   solsize = cell2mat(cdinfo.nofgrid);
   p = sum(solsize(1:3))+24;


   [ia,ja,nzA] = find(A);
   nA = sparse(ja,1:length(ja),1); 


   c1            = receive(socketp);
   c2            = receive(socketp);
   c3            = receive(socketp);
   c4            = receive(socketp);
   c5            = receive(socketp);
   c6            = receive(socketp);

   c1r = nA*(c1.*nzA.*dlambda); 
   c2r = nA*(c2.*nzA.*dlambda); 
   c3r = nA*(c3.*nzA.*dlambda); 
   c4r = nA*(c4.*nzA.*dlambda); 
   c5r = nA*(c5.*nzA.*dlambda); 
   c6r = nA*(c6.*nzA.*dlambda);

   send(socketp,c1r);
   send(socketp,c2r);
   send(socketp,c3r);
   send(socketp,c4r);
   send(socketp,c5r);
   send(socketp,c6r);

   dA             = receive(socketp);
  
   yf            = receive(socketp);
   zf            = receive(socketp);  
   y10f            = receive(socketp);
   z10f            = receive(socketp);
   y01f            = receive(socketp);
   z01f            = receive(socketp);

   
 
%figure
%hold on
%plot(yf,zf,'or');
%plot(y10f,z10f,'xg')
%plot(y01f,z01f,'xb')

wlist = find(sum(A')<0.9);

%h = plot(yf(wlist),zf(wlist),'or');
%set(h,'linewidth',3)



   
   dA = dA(1:end-24)';
   dA = [dA ,zeros(1,cdinfo.nofgrid{4})];

   if(min(sum(A'))<0.5)
      disp('A is wrong !!!!!!!')
   end


 
   closeport(socketp);
