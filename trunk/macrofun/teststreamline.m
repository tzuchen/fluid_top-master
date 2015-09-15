% This is a routine to test the streamline computation
%
% Tzu-Chen Liang   11/16/2005

% x has to be precomputated by solver 
% and also we need cdinfo
% sol = x;
ds  = 10;

[Y,Z] = meshgrid(0.1:0.1:0.9,0.1:0.1:0.9);
[m,n] = size(Y);

 x0   = [0.01 *ones(1,m*n); 
          reshape(Y,1,m*n) ;
          reshape(Z,1,m*n)];  

 x       = x0;
 xstream =[];
 ncount  = 1;
 xstream(:,:,1)= x0;
 ncount = 0;
    while and(min(x(1,:))<cdinfo.coord{1}.r,ncount<500)
       ncount=ncount + 1;
       [x,w] = streamlinestep(x,sol,ds,cdinfo);
       xstream(:,:,ncount) = x;
    end

 for ncount = 1:m*n
      sline = squeeze(xstream(:,ncount,:));
      plot3(sline(1,:),sline(2,:),sline(3,:));
      hold on
 end

 x = x(2:3,:);
 Am = lininterp2dmat(x,Y,Z);
