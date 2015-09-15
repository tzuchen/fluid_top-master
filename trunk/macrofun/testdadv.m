% Test routine to find dA/dv 
%
% Tzu-Chen Liang 2-13-2006

 
        indu = 1: cdinfo.nofgrid{1}; 
        indv = cdinfo.nofgrid{1}+1: cdinfo.nofgrid{1}+cdinfo.nofgrid{2}; 
        indw = cdinfo.nofgrid{1}+cdinfo.nofgrid{2}+1: cdinfo.nofgrid{1}+cdinfo.nofgrid{2}+cdinfo.nofgrid{3}; 
 
        vellen   = cdinfo.nofgrid{1}+ cdinfo.nofgrid{2}+ cdinfo.nofgrid{3};         

        
 
 

   	% generate a set of point to begin the streamlines.
   	dy = 0.1/2;
        yspam = 0:dy:1;
        zspam = 0:dy:1;
        ny    = length(yspam);
        nz    = length(zspam);
   	[Y0,Z0] = meshgrid(yspam,zspam);
        
   	nofline = prod(size(Y0));
   	x0 = (cdinfo.coord{1}.l)*ones(1,nofline);
   	y0 = reshape(Y0,1,nofline); 
   	z0 = reshape(Z0,1,nofline);
   	y  = y0;
   	z  = z0;
        n  = length(y);
 
        vlen = 0.01;
        y10 = y + vlen;
        z10 = z;
        y01 = y;
        z01 = z + vlen;
  
        yf   = y;
        zf   = z;
        y10f = y10;
        z10f = z10;
        y01f = y01; 
        z01f = z01;


     
        [yf  ,zf  ,Df  ] = streamlinemap(yf  ,zf  ,1,cdinfo,sol);
        [y10f,z10f,D10f] = streamlinemap(y10f,z10f,1,cdinfo,sol);
        [y01f,z01f,D01f] = streamlinemap(y01f,z01f,1,cdinfo,sol);
   
        Dyf   = Df.dydv;
        Dzf   = Df.dzdv;
        Dy10f = D10f.dydv;
        Dz10f = D10f.dzdv;
        Dy01f = D01f.dydv;
        Dz01f = D01f.dzdv;

 
        p11  = -Dzf + Dz01f;
        p12  =  Dzf + Dz10f;
        p21  =  Dyf - Dy01f;
        p22  = -Dyf + Dy10f; 
           
        p1 = (y10f - yf);
        q1 = (z10f - zf);
        p2 = (y01f - yf);
        q2 = (z01f - zf);

        A  = sparse(n,n);
        dA = sparse(n^2,vellen);
        
        for i = 1:n
                 
           Pmat   = [p1(i) p2(i); q1(i) q2(i)]/vlen;
           invP   = inv(Pmat);
 
           pos    = invP*[y-yf(i);z-zf(i)];
    

           DinvP  = 1/det(Pmat)*[p11(i,:);p21(i,:);p12(i,:);p22(i,:)];

           i
 

           for j = 1:n
               d = norm(pos(:,j));
               if abs(d)<=2*dy
                      A(j,i) = exp(-(d/dy).^2)/pi;
                      ij     = sub2ind(size(A),i,j);
                      coeff  = reshape(pos(:,j)*[y(j)-yf(i)  z(j)-zf(i)],1,4);
                      
                      dA(ij,:) = A(j,i)/pi*2/(dy^2)*(coeff*DinvP - pos(:,j)'*invP*[Dyf(i,:);Dzf(i,:)]);
                      
        


 
               end
           end
        end

        for i = 1:n
                   A(i,:) = A(i,:)/sum(A(i,:));
        end
 







 

	%X1 = zeros(ny,nz);
	%X1((ny-1)/2+1:end,:) = 1;
	%x1 = reshape(X1,(ny)*(nz),1);
 
