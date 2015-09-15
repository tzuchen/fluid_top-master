%
% Given a map, this routine evolves a line
% We can conclude the length of the line grows linearly
% or exponentionally
% 
%
%



mapfunction = @standardmap
param{1} = 0.006;
period   = 1;
 
%mapfunction = @twistmap;
%param{1}    = 0.5;
%param{2}    = 0.5;
%param{3}    = 1;
%period      = 0;

%mapfunction = @linkedtwistmap;
%param{1}    = 0.5;
%param{2}    = 0.4;
%param{3}    = 1;
%param{4}    = 0.5
%param{5}    = 0.6
%param{6}    = 1
%period      = 0;


nic          = 60;
showmovie    = 0;

close all


bclen = [];
% generate points on a line to begin the streamlines.
  dy      = 0.02;
  [y0,z0] = meshgrid(0:dy:1,0.5);
  nofline = length(y0);
  y       = y0;
  z       = z0;
  
  %blen    = sum(sqrt(diff(y).^2 + diff(z).^2));
  blen    = sum(sqrt(min(diff(y),1-diff(y)).^2 + min(diff(z),1-diff(y)).^2));
       
  clear M
 if showmovie ==1
  %figure
  %hold on
 end
  
  for i= 1:nic

      dsq    = sqrt(diff(y).^2 + diff(z).^2);
      t      = [0 cumsum(dsq)];
      bratio = t(end)/blen;
             
      pt     = ceil(nofline*bratio);
      y      = interp1(t,y,(0:pt)/pt*t(end));
      z      = interp1(t,z,(0:pt)/pt*t(end));
      
      if i>1 
           %y    = y(2:end-1);
           %z    = z(2:end-1);   
  
           [y,z] = feval(mapfunction,y,z,1,1,param);  
           
              
           %y     = [y0(1) y y0(end)];
           %z     = [z0(1) z z0(end)];
       end   
      
      disp(sprintf('count : %d  length: %f',i,t(end)))
      bclen = [bclen t(end)];
      
      if showmovie ==1       
         %fill([0 y 1],[0 z 0],'r');
         %fill([0 y 1],[1 z 1],'b');  
          plot(y,z,'.')             
         axis equal
         axis([0 1 0 1])     
         M(i) = getframe;
      end
  end
      
      if showmovie ==1
         %fill([0 y 1],[0 z 0],'r');
         %fill([0 y 1],[1 z 1],'b');
         axis equal
      end

figure
plot(bclen);

