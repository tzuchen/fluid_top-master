% This is a demo of the function streamlinemap.m 
%
% Notice:
%   make sure sol and cdinfo exist
%
% Tzu-Chen Liang  1/12/2006
   close all

demo = 2;  

% make the pressure part zero
sol(end-cdinfo.nofgrid{4}+1:end)=0;


switch demo

%
% Demo 0: just calculate the streamline and the mapping
%
case 0
   	% generate a set of point to begin the streamlines.
   	dy = 0.1;
   	[Y0,Z0] = meshgrid(0.1:dy:0.9,0.1:dy:0.9);
   	nofline = prod(size(Y0));
   	x0 = (cdinfo.coord{1}.l)*ones(1,nofline);
   	y0 = reshape(Y0,1,nofline); 
   	z0 = reshape(Z0,1,nofline);
   	y  = y0;
   	z  = z0;

   	
        [y,z,S,sm] = streamlinemap(y0,z0,cdinfo,sol);
 
        figure 
        hold on             
          
        for i=1:size(sm.x,2)
          plot3(sm.x(:,i),sm.y(:,i),sm.z(:,i))
        end
try, gdplot(alinp',4,cdinfo,[param(5) color(5)]),end

%
% Demo 1: for regular grids, calculate the mapping, repeat n times.
%
case 1
   	% generate a set of point to begin the streamlines.
   	dy = 0.1/2;
   	[Y0,Z0] = meshgrid(0:dy:1,0:dy:1);
   	nofline = prod(size(Y0));
   	x0 = (cdinfo.coord{1}.l)*ones(1,nofline);
   	y0 = reshape(Y0,1,nofline); 
   	z0 = reshape(Z0,1,nofline);
   	y  = y0;
   	z  = z0;

   	figure
   	hold on
   	h = plot(y,z,'o');
   	set(h,'linewidth',2)
  
   	for i= 1:1
              [y,z,spm] = streamlinemap(y,z,cdinfo,sol);
   	      h = plot(y,z,'xr');
   	      set(h,'linewidth',2);
   	end
        for i=1:nofline
            h= line([y0(i),y(i)],[z0(i),z(i)]);
            set(h,'linewidth',2,'color','g');
        end
   	axis equal
        axis([0 1 0 1])
%
% Demo 2: Evolve the boundry between red/blue liquid, repeat n times   
%
case 2
       % generate points on a line to begin the streamlines.
       dy = 0.02;
       [y0,z0] = meshgrid(0:dy:1,0.5);
       nofline = length(y0);
       x0 = (cdinfo.coord{1}.l)*ones(1,nofline);
       y  = y0;
       z  = z0;
       blen = sum(sqrt(diff(y).^2 + diff(z).^2));
       clear M
   	figure
   	hold on
   	%h = plot(y,z,'o');
   	%set(h,'linewidth',2)
  
 
   	for i= 1:300


              dsq   = sqrt(diff(y).^2 + diff(z).^2);
              t     = [0 cumsum(dsq)];

              bratio = t(end)/blen;
             
              pt = ceil(nofline*bratio);
              y = interp1(t,y,(0:pt)/pt*t(end));
   	      z = interp1(t,z,(0:pt)/pt*t(end));
          if i>1 
              y = y(2:end-1);
              z = z(2:end-1);   
     
                           
                 [y,z,spm] = streamlinemap(y,z,cdinfo,sol);
              
              y     = [0; y; 1]';
              z     = [0.5; z; 0.5]';
          end


              %h = plot(y,z,'xr');
              %h = plot(y,z,'g');    
     
      
              disp(sprintf('count : %d',i))
              
              fill([0 y 1],[0 z 0],'r');
              fill([0 y 1],[1 z 1],'b');  
             
              axis equal
              axis([0 1 0 1])     
              M(i) = getframe;
   	end
              fill([0 y 1],[0 z 0],'r');
              fill([0 y 1],[1 z 1],'b');
              save M
              %h = plot(y,z,'xr');
   	      %set(h,'linewidth',2)
              %h = plot(y,z,'r'); 
              %set(h,'linewidth',2)     
        axis equal
%
% Demo 3: calculate the stretch of the grids 
%
case 3
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

   	figure
   	hold on
   	h = plot(y,z,'o');
   	set(h,'linewidth',2)
        yf = y;
        zf = z;
        y10f = y10;
        z10f = z10;
        y01f = y01; 
        z01f = z01;
     
  
           for i = 1:1
              [yf  ,zf  ,spmf  ] = streamlinemap(yf  ,zf  ,cdinfo,sol);
              [y10f,z10f,spm10f] = streamlinemap(y10f,z10f,cdinfo,sol);
              [y01f,z01f,spm01f] = streamlinemap(y01f,z01f,cdinfo,sol);
           end
 
   	      h = plot(yf,zf,'xr');
   	      set(h,'linewidth',2)
              h = plot(y10f,z10f,'xr');
              h = plot(y01f,z01f,'xr');

             for i = 1:length(y)
                   plot([yf(i),y10f(i)],[zf(i),z10f(i)]);
                   plot([yf(i),y01f(i)],[zf(i),z01f(i)]);
             end   	
       
               p1 = (y10f - yf);
               q1 = (z10f - zf);
               p2 = (y01f - yf);
               q2 = (z01f - zf);

               A  = sparse(n,n);
               dA = sparse(n^2,size(spmf,2));
               for i = 1:n
                 
                 Pmat = [p1(i) p2(i); q1(i) q2(i)]/vlen;
                 pos  = inv(Pmat)*[y-yf(i);z-zf(i)];


                 for j = 1:n
                   d = norm(pos(:,j));
                   if abs(d)<=2*dy
                      A(j,i) =exp(-(d/dy).^2)/pi;


                   end 
		 end
               end
              for i = 1:n
                   A(i,:) = A(i,:)/sum(A(i,:));
              end
 
   	axis equal
 

	X1 = zeros(ny,nz);
	X1((ny-1)/2+1:end,:) = 1;
	x1 = reshape(X1,(ny)*(nz),1);
        figure
	contourf(reshape((A^1)*x1,ny,nz))
	figure
	contourf(reshape((A^10)*x1,ny,nz))
	%figure
	%contourf(reshape((A^100)*x1,ny,nz)) 


%
% Demo 4 : Calculate the linear approximation of the finite resolution mapping 
%
case 4
   	% generate a set of point to begin the streamlines.
   	dy = 1/50;
        yspan = 0:dy:1;
        zspan = 0:dy:1;

        ny    = length(yspan);
        nz    = length(zspan);
   	[Y0,Z0] = meshgrid(yspan,zspan);
   	nofline = prod(size(Y0));
   	x0 = (cdinfo.coord{1}.l)*ones(1,nofline);
   	y0 = reshape(Y0,1,nofline); 
   	z0 = reshape(Z0,1,nofline);
   	y  = y0;
   	z  = z0;

   	figure
   	hold on
   	h = plot(y,z,'o');
   	set(h,'linewidth',2)
  
   	for i= 1:1
              [y,z,spm] = streamlinemap(y,z,cdinfo,sol);
   	      h = plot(y,z,'xr');
   	      set(h,'linewidth',2)
   	end

       
        Y = reshape(y,ny,nz);
        Z = reshape(z,ny,nz);
        
   	axis equal

case 5
   	dy = 0.05;
        yspan = 0.5;
        zspan = 0.3:dy:0.97;

        ny    = length(yspan);
        nz    = length(zspan);
   	[Y0,Z0] = meshgrid(yspan,zspan);
   	nofline = prod(size(Y0));
   	x0 = (cdinfo.coord{1}.l)*ones(1,nofline);
   	y0 = reshape(Y0,1,nofline); 
   	z0 = reshape(Z0,1,nofline);

        %y0(find(mod(1:nofline,100)~=1))=0;
        %z0(find(mod(1:nofline,100)~=1))=0;

   	y  = y0;
   	z  = z0;

   	figure
   	hold on
   	h = plot(y,z,'o');
   	set(h,'linewidth',2)
  
   	for i= 1:1
              [y,z,spm] = streamlinemap(y,z,cdinfo,sol);
   	      h = plot(y,z,'xr');
   	      set(h,'linewidth',2)
   	end

       
      %  Y = reshape(y,ny,nz);
      %  Z = reshape(z,ny,nz);
        
   	axis equal

%
% Demo 6: Examine that if the mix-norm dorps a fixed factor after mapping   
%
case 6
       % generate points on a line to begin the streamlines.
       dy = 0.02;
       [y0,z0] = meshgrid(0:dy:1,0.5);
       nofline = length(y0);
       x0 = (cdinfo.coord{1}.l)*ones(1,nofline);
       y  = y0;
       z  = z0;
       blen = sum(sqrt(diff(y).^2 + diff(z).^2));
   
   	figure
   	hold on
   	h = plot(y,z,'o');
   	set(h,'linewidth',2)
  
 
   	for i= 1:600


              dsq   = sqrt(diff(y).^2 + diff(z).^2);
              t     = [0 cumsum(dsq)];

              bratio = t(end)/blen;
             
              pt = ceil(nofline*bratio);
              y = interp1(t,y,(0:pt)/pt*t(end));
   	      z = interp1(t,z,(0:pt)/pt*t(end));

              y = y(2:end-1);
              z = z(2:end-1);   
              
              [y,z,spm] = streamlinemap(y,z,cdinfo,sol);
              y     = [0, y, 1];
              z     = [0.5, z, 0.5];
              yb{i} = y;
              zb{i} = z;

   	end
       ind = 0;     
       sigma=[];
       for i = 1:10:600
          ind = ind+1;
          sigma(ind) = mixnormI(yb{i},zb{i},20,1)
       end

end














