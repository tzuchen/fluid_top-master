clear all

nproc       =  5;
den         =  5;
bodyf{1}    =  1;
bodyf{2}    =  0;
bodyf{3}    =  0;
mu          =  1;
ny          = 50;
nz          = ny;

[Ap,Ac,Pp,Pr,bp,p2Amatp,p,r,cdinfo] = problemgen(nproc,den,bodyf,mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set alpha
  clear alphapos
%                xl   xu   yl   yu   zl   zu  
  alphapos{1} = [0.5   1    0.2    0.5    0.2   0.5];
  alphapos{2} = [0.5   1    0.5    0.8    0.5   0.8]; 
  alphapos{3} = [1    1.5   0.2    0.5    0.5   0.8];  
  alphapos{4} = [1    1.5   0.5    0.8    0.2   0.5];
  alphaind    = alphaset(alphapos,cdinfo);
  al2pmat     = al2p(alphaind,cdinfo);


  %alphaijk  = ind2ijk(1:cdinfo.nofgrid{4},4,cdinfo);
  %alphaind{1} = find(and(abs(alphaijk(1,:)-10-alphaijk(2,:))<2,alphaijk(3,:)<=4));  
  %al2pmat  = al2p(alphaind,cdinfo); 
  
  alpha       = 10000*ones(length(cell2mat(alphaind)),1);
  alphainp    = al2pmat*alpha;
  alphavec    = p2Amatp*alphainp; 
   
  cp = bp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydim    = cdinfo.coord{2}.r - cdinfo.coord{2}.l;
zdim    = cdinfo.coord{3}.r - cdinfo.coord{3}.l;
dy      = ydim / ny;
dz      = zdim / nz;
ypt     = dy/2:dy:ydim-dy/2;
zpt     = dz/2:dz:zdim-dz/2;

[Y0,Z0] = meshgrid(ypt,zpt);
nofline = prod(size(Y0)); 

x0      = (cdinfo.coord{1}.l)*ones(1,nofline);
y0      = reshape(Y0,1,nofline); 
z0      = reshape(Z0,1,nofline);

xcor    = (cdinfo.coord{1}.l)*ones(1,(ny+1)*(nz+1));
[Yc,Zc] = meshgrid(0:dy:ydim,0:dz:zdim);
ycor    = reshape(Yc,(ny+1)*(nz+1),1);
zcor    = reshape(Zc,(ny+1)*(nz+1),1);




solsize =  sum(cell2mat(cdinfo.nofgrid));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on 

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now the petsc part
  PetscInitialize(nproc,'taskmanager');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    socketp = openport();
 
    send(socketp,Ap);
    send(socketp,Ac);
    send(socketp,Pp);
    send(socketp,bp);
    send(socketp,cp);
    
    send(socketp,x0); 
    send(socketp,y0); 
    send(socketp,z0); 

    send(socketp,xcor); 
    send(socketp,ycor); 
    send(socketp,zcor); 

    send(socketp,dy);
    send(socketp,dz);
  
    cdinfosend(socketp,cdinfo);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %loop here
    for ncount = 1:5
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Task 1: solve the velocity field
           disp('Task 1')
          send(socketp,1); 

   	    send(socketp,diag(sparse(alphavec))); 
   	    xp    = receive(socketp); % Receive xp
            x     = xp(r); 
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Task 2: Solve A
          disp('Task 2')
          send(socketp,2); 

            sol  = x; 
            sol(end-cdinfo.nofgrid{4}+1:end)=0;
            send(socketp,sol);
            send(socketp,ny);
            send(socketp,nz);  


            OPTS.disp = 0;

            % just for it to function correctly
            [A]         = SparseReceive(socketp,ny*nz,ny*nz); 
            A           = fixA(A);
            [V,D]       = eigs(A',6,'LM',OPTS);
            [Dr,I]      = sort(max(abs(D)));
            Im          = I(end-2);
            mu          = D(I(end-2),I(end-2));
            nzInd       = find(A);
            dlambda     = real(V(:,Im))*real(V(:,Im))';
            dlambda     = dlambda(nzInd);

       
           


            [ia,ja,nzA] = find(A);
            nA = sparse(ja,1:length(ja),1); 


            c1            = receive(socketp);
   	    c2            = receive(socketp);
   	    c3            = receive(socketp);
   	    c4            = receive(socketp);
   	    c5            = receive(socketp);
   	    c6            = receive(socketp);


          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Task 3: solve dA
          disp('Task 3')
          send(socketp,3); 

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
          dA             = dA(1:end-24)';
          dA             = [dA ,zeros(1,cdinfo.nofgrid{4})];
          cp             = dA(p);

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Task 5: getA and dA
          disp('Task 5')
          send(socketp,5);

            sol  = x; 
            sol(end-cdinfo.nofgrid{4}+1:end)=0;
            send(socketp,sol);
             
            y = receive(socketp);
            z = receive(socketp);
	    Da.dydv = SparseReceive(socketp,ny*nz,solsize);
	    Da.dzdv = SparseReceive(socketp,ny*nz,solsize);
	    yecor = receive(socketp);
	    zecor = receive(socketp);


            Ye    = reshape(y    ,ny  ,nz);
            Ze    = reshape(z    ,ny  ,nz);
            Ycor  = reshape(yecor,ny+1,nz+1);
            Zcor  = reshape(zecor,ny+1,nz+1);
          



          [A,dA]=AdAfind(Y0,Z0,Ye,Ze,Ycor,Zcor,ny,nz,sol,cdinfo,Da);
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Task 4: solve gradient
          disp('Task 4')
          send(socketp,4);

            send(socketp,diag(sparse(alphavec)));
            send(socketp,cp);
            send(socketp,xp);

   	    fobj  = receive(socketp); % Receive fobj
   	    fobj  = fobj(1);  
   	    gradf = receive(socketp); % Receive gradf

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

           gradf = gradf/max(abs(gradf));                            
           alphainp = alphainp + gradf*2000;          
           alphavec = p2Amatp*alphainp;
          %alphavec = alphavec  + p2Amatp*gradf*1000;
           alphavec(find(alphavec<0)) = 0;

           mu
         %  h = plot(sort(diag(abs(D))),'r');
         %  set(h,'linewidth',2);
    end

    
   PetscFinalize(socketp);
 





sol = x;
% we don't need pressure
sol(end-cdinfo.nofgrid{4}+1:end)=0;













