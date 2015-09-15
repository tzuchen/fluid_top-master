% Dem to evolve the boundry
%
% Tzu-Chen Liang  3-29-2006

iter = 5;

% generate points on a line to begin the streamlines.
  
  dy = 0.02;
  [y0,z0] = meshgrid(0:dy:1,0.5);
  nofline = length(y0);
  x0 = (cdinfo.coord{1}.l)*ones(1,nofline);
  y  = y0;
  z  = z0;
  blen = sum(sqrt(diff(y).^2 + diff(z).^2));
  clear Mb
  figure
  hold on

  for i= 1:iter

              dsq   = sqrt(diff(y).^2 + diff(z).^2);
              t     = [0 cumsum(dsq)];

              bratio = t(end)/blen;
             
              pt = ceil(nofline*bratio);
              y = interp1(t,y,(0:pt)/pt*t(end));
   	      z = interp1(t,z,(0:pt)/pt*t(end));

              y = y(2:end-1);
              z = z(2:end-1);   
              
              [y,z,spm] = streamlinemap(y,z,cdinfo,sol);
              y     = [0; y; 1]';
              z     = [0.5; z; 0.5]';   
      
              disp(sprintf('count : %d',i))
              
              fill([0 y 1],[0 z 0],'r');
              fill([0 y 1],[1 z 1],'b');  
             
              axis equal
              axis([0 1 0 1])     
              Mb(i) = getframe;
 end
 
 %save Mb
