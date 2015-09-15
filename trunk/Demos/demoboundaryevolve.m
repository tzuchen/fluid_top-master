

load YZ400

[Y0,Z0]=meshgrid(0:1/400:1,0:1/400:1);

y = [0:0.1:1];
z = 0*y+0.5;

iter=150


  for i= 1:iter

              dsq   = sqrt(diff(y).^2 + diff(z).^2);
              t     = [0 cumsum(dsq)];

              bratio = t(end)/blen;
             
              pt = ceil(nofline*bratio);
              y = interp1(t,y,(0:pt)/pt*t(end));
   	      z = interp1(t,z,(0:pt)/pt*t(end));

              y = y(2:end-1);
              z = z(2:end-1);   
              
              y = interp2(Y0,Z0,Ye,y,z);
              z = interp2(Y0,Z0,Ze,y,z);
              y     = [0 y 1];
              z     = [0.5 z 0.5];   
      
              disp(sprintf('count : %d',i))
              
              fill([0 y 1],[0 z 0],'r');
              fill([0 y 1],[1 z 1],'b');  
             
              axis equal
              axis([0 1 0 1])     
              Mb(i) = getframe;
 end
