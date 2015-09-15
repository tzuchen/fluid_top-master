% democutoffdef
close all
ptsize = 20;
h =[];
k=1;
d = 0.01;
nc =4
 for p = 0:4
  for q = 0:3
      
      
        %h(p+1,q+1) = axes('Position',[ 0.08+p/5 0.90*(q/nc+1*d)  0.9*(1/5) 0.90*(1/nc-d) ]);
        h(p+1,q+1) = axes('Position',[ 0.08+0.9*(p/5) 0.05+0.95*(q/nc)  0.9*(1/5) 0.95*(1/nc) ]);
        axis equal
        
        k = k+1;  
  end
 end
 
 axisvec = [-0.2 1.2 -0.2 1.2 -0.2 1.2]
 
 for i = 1:4
     axes(h(4,i))
     if i==4 
        fr =scatter([0 1 2],[-0.4 -0.4 -0.4],20,[0 0 0],'filled')
     else
       if i==3
        fr =scatter([0 1 2],[-0.2 -0.2 -0.2],20,[0 0 0],'filled')
       else
       fr =scatter([0 1 2],[0 0 0],20,[0 0 0],'filled')     
       end
     end
       % set(fr,'color',[0 0 0])
       axis([-0.5 2.5 -0.5 0.5]) 
       grid off
       box off
              set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'visible','off')
 end
 
 h = h(:,end:-1:1)
 
% 0d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i = 1:5
%    if i~=4
%       axes(h(i,1))
%       X = [0];
%       Y = [0];
%       Z = [0]
%       w =  scatter3(X,Y,Z,ptsize,'filled');
%       %set(w,'markersize',5,'linewidth',2) 
%       grid off
%        set(gca,'xtick',[]);
%        set(gca,'ytick',[]);
%        set(gca,'ztick',[]);
%        box off
%        axis(axisvec) 
%        hold on
%        %plot3([0 0],[0 0],[0 1],'r','linewidth',2)
%        w =  scatter3(X,Y,Z,400,'filled','r');
%        set(gca,'visible','off')
%    end   
%end 
 
% 1 d
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for i = 1:5
    dz = -0.3
    if i~=4
       axes(h(i,1))
       X = [0 0];
       Y = [0 1];
       Z = [0+dz 0+dz];
       if i==1
       w =  scatter3(X,Y,Z,ptsize,'filled');
       end
       hold on 
       %set(w,'markersize',5,'linewidth',2) 
       grid off
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'ztick',[]);
        box off
        axis(axisvec) 
        plot3(X,Y,Z)
        
        if i==1
           %plot3([0 0],[0 0],[0 1],'r','linewidth',2) 
            w =  scatter3(0,0,0+dz,400,'filled','r');
        end
 
        
          if i >1
            for j = 1:2
                % plot3([X(j) X(j)],[Y(j) Y(j)],[Z(j) Z(j)+0.5],'r','linewidth',2);
                 w =  scatter3(X,Y,Z,200,'filled','r');           
            end
          end        
        set(gca,'visible','off')
       %gh = plot3(X,Y,Z)
      % set(gh,'linewidth',2);
        set(gca,'visible','off')
        view(3)
    end   
end
%2d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i = 1:5
   
    if i~=4
         axes(h(i,2))
       X = [0 0 1 1]; 
       Y = [0 1 1 0];
       Z = [0 0 0 0];
       
       hold on 
       %set(w,'markersize',5,'linewidth',2) 
       grid off
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'ztick',[]);
        box off
        axis(axisvec) 
         plot3([0 0],[0 1],[0 0])
         plot3([0 1],[1 1],[0 0])
         plot3([1 1],[1 0],[0 0])
         plot3([1 0],[0 0],[0 0])
        
        if i==1

           %plot3([0 0],[0 0],[0 1],'r','linewidth',2) ;
            w =  scatter3(0,1,0,ptsize,'filled','b');
            w =  scatter3(1,0,0,ptsize,'filled','b');
            w =  scatter3(1,1,0,ptsize,'filled','b');
            w =  scatter3(0,0,0,400,'filled','r'); 
        end
        if i ==2
            w =  scatter3(1,1,0,ptsize,'filled','b');
            %  plot3([0 0],[0 0],[0 1/3],'r','linewidth',2) ;
            %  plot3([0 0],[1 1],[0 1/3],'r','linewidth',2) ;
            %  plot3([1 1],[0 0],[0 1/3],'r','linewidth',2) ;
            w =  scatter3(0,0,0,133,'filled','r'); 
            w =  scatter3(0,1,0,133,'filled','r');
            w =  scatter3(1,0,0,133,'filled','r');
                          
        end    
        if i ==3
              %plot3([0 0],[0 0],[0 8/27],'r','linewidth',2) ;
              %plot3([0 0],[1 1],[0 19/81],'r','linewidth',2) ;
              %plot3([1 1],[0 0],[0 19/81],'r','linewidth',2) ;
              %plot3([1 1],[1 1],[0 19/81],'r','linewidth',2) ;     
              w =  scatter3(0,0,0,120,'filled','r'); 
              w =  scatter3(0,1,0,90,'filled','r');
              w =  scatter3(1,0,0,90,'filled','r');
              w =  scatter3(1,1,0,40,'filled','r');
              
        end 
        if i ==5
            for j = 1:4
                 %plot3([X(j) X(j)],[Y(j) Y(j)],[Z(j) Z(j)+0.25],'r','linewidth',2);
                 w =  scatter3(X,Y,Z,100,'filled','r');
            end
        end   
        view(3)
        set(gca,'visible','off')
    end   
 end
 
 
%3d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i = 1:5
   
    if i~=4
         axes(h(i,3))
       X = [0 0 1 1 0 0 1 1]; 
       Y = [0 1 1 0 0 1 1 0];
       Z = [0 0 0 0 1 1 1 1];
       w =  scatter3(X,Y,Z,ptsize,'filled');
       hold on 
       plot3([0 0],[0 1],[0 0])
       plot3([0 1],[1 1],[0 0])
       plot3([1 1],[1 0],[0 0])
       plot3([1 0],[0 0],[0 0])
       plot3([0 0],[0 1],[1 1])
       plot3([0 1],[1 1],[1 1])
       plot3([1 1],[1 0],[1 1])
       plot3([1 0],[0 0],[1 1])
       plot3([0 0],[0 0],[0 1])
       plot3([1 1],[1 1],[0 1])
       plot3([1 1],[0 0],[0 1])
       plot3([0 0],[1 1],[0 1])
       
      % set(w,'markersize',5,'linewidth',2) 
       grid off
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'ztick',[]);
        box off
        axis(axisvec) 
        if i==1
          % plot3([0 0],[0 0],[1 2],'r','linewidth',2) 
            w =  scatter3(0,0,0,400,'filled','r'); 
        end
        if i ==2
              %plot3([0 0],[0 0],[1 1+1/4],'r','linewidth',2) ;
              %plot3([1 1],[0 0],[1 1+1/4],'r','linewidth',2) ;
              %plot3([0 0],[1 1],[1 1+1/4],'r','linewidth',2) ;
              %plot3([0 0],[0 0],[0 1/4],'r','linewidth',2) ;  
              
              w =  scatter3(0,0,0,140,'filled','r');
              w =  scatter3(0,1,0,140,'filled','r');
              w =  scatter3(0,0,1,140,'filled','r');
              w =  scatter3(1,0,0,140,'filled','r');
              
        end    
        if i ==3
              %plot3([0 0],[0 0],[0 8/27],'r','linewidth',2) ;
              %plot3([0 0],[1 1],[0 19/81],'r','linewidth',2) ;
              %plot3([1 1],[0 0],[1 1+19/81],'r','linewidth',2) ;
              %plot3([1 1],[1 1],[1 1+19/81],'r','linewidth',2) ;   
              %plot3([1 1],[0 0],[0 1/4],'r','linewidth',2) ;
              %plot3([0 0],[0 0],[1 1+8/27],'r','linewidth',2) ;
              %plot3([0 0],[1 1],[1 1+19/81],'r','linewidth',2) ;
              %plot3([1 1],[0 0],[1 1+19/81],'r','linewidth',2) ;
              
              w =  scatter3(0,0,0,120,'filled','r');
              w =  scatter3(0,1,0,100,'filled','r');
              w =  scatter3(0,0,1,100,'filled','r');
              w =  scatter3(1,0,0,100,'filled','r');
              w =  scatter3(1,1,0,60,'filled','r');
              w =  scatter3(0,1,1,60,'filled','r');
              w =  scatter3(1,0,1,60,'filled','r');
              %w =  scatter3(1,0,0,25,'filled','r');
              
        end 
        if i ==5
            for j = 1:8
                 %plot3([X(j) X(j)],[Y(j) Y(j)],[Z(j) Z(j)+0.125],'r','linewidth',2);
                 w =  scatter3(X,Y,Z,50,'filled','r');
            end
        end   
        set(gca,'visible','off')
    end   
  end
 %4d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for i = 1:5
     dx = 0.13
     dy = 0.13
     dz = 0.13
    if i~=4
        axes(h(i,4))
       X = [0 0 1 1 0 0 1 1]; 
       Y = [0 1 1 0 0 1 1 0];
       Z = [0 0 0 0 1 1 1 1];
       X = [X-dx X+dx];
       Y = [Y-dy Y+dy];
       Z = [Z-dz Z+dz];
       
       

       
        w =  scatter3(X,Y,Z,ptsize,'filled');
       hold on 
       plot3([0 0]+dx,[0 1]+dy,[0 0]+dz)
       plot3([0 1]+dx,[1 1]+dy,[0 0]+dz)
       plot3([1 1]+dx,[1 0]+dy,[0 0]+dz)
       plot3([1 0]+dx,[0 0]+dy,[0 0]+dz)
       plot3([0 0]+dx,[0 1]+dy,[1 1]+dz)
       plot3([0 1]+dx,[1 1]+dy,[1 1]+dz)
       plot3([1 1]+dx,[1 0]+dy,[1 1]+dz)
       plot3([1 0]+dx,[0 0]+dy,[1 1]+dz)
       plot3([0 0]+dx,[0 0]+dy,[0 1]+dz)
       plot3([1 1]+dx,[1 1]+dy,[0 1]+dz)
       plot3([1 1]+dx,[0 0]+dy,[0 1]+dz)
       plot3([0 0]+dx,[1 1]+dy,[0 1]+dz)
       
       plot3([0 0]-dx,[0 1]-dy,[0 0]-dz)
       plot3([0 1]-dx,[1 1]-dy,[0 0]-dz)
       plot3([1 1]-dx,[1 0]-dy,[0 0]-dz)
       plot3([1 0]-dx,[0 0]-dy,[0 0]-dz)
       plot3([0 0]-dx,[0 1]-dy,[1 1]-dz)
       plot3([0 1]-dx,[1 1]-dy,[1 1]-dz)
       plot3([1 1]-dx,[1 0]-dy,[1 1]-dz)
       plot3([1 0]-dx,[0 0]-dy,[1 1]-dz)
       plot3([0 0]-dx,[0 0]-dy,[0 1]-dz)
       plot3([1 1]-dx,[1 1]-dy,[0 1]-dz)
       plot3([1 1]-dx,[0 0]-dy,[0 1]-dz)
       plot3([0 0]-dx,[1 1]-dy,[0 1]-dz)
       
       plot3([0-dx 0+dx],[0-dy 0+dy],[0-dz 0+dz])
       plot3([0-dx 0+dx],[0-dy 0+dy],[1-dz 1+dz])
       plot3([0-dx 0+dx],[1-dy 1+dy],[0-dz 0+dz])
       plot3([1-dx 1+dx],[0-dy 0+dy],[0-dz 0+dz])
       plot3([1-dx 1+dx],[1-dy 1+dy],[0-dz 0+dz])
       plot3([0-dx 0+dx],[1-dy 1+dy],[1-dz 1+dz])
       plot3([1-dx 1+dx],[0-dy 0+dy],[1-dz 1+dz])
       plot3([1-dx 1+dx],[1-dy 1+dy],[1-dz 1+dz])
             
       
       
       %set(w1,'markersize',5,'linewidth',2)
       %set(w2,'markersize',5,'linewidth',2)
       
       grid off
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'ztick',[]);
        box off
        axis(axisvec) 
        if i==1
           %plot3([0 0],[0 0],[1 2],'r','linewidth',2) 
           w =  scatter3(0-dx,0-dx,0-dx,400,'filled','r');
           text(0.2,-0.3,-0.3,'$\omega_4^0$','fontsize',16,'interpreter','latex')
        end
        if i ==2
              %plot3([0 0],[0 0],[1 1+1/5],'r','linewidth',2) ;
              %plot3([1 1],[0 0],[1 1+1/5],'r','linewidth',2) ;
              %plot3([0 0],[1 1],[1 1+1/5],'r','linewidth',2) ;
              %plot3([0 0],[0 0],[0 1/5],'r','linewidth',2) ;  
              %plot3([0+0.4 0+0.4],[0 0],[0.3 0.3+1/5],'r','linewidth',2) ; 
              
              w =  scatter3(0-dx,0-dy,1-dz,100,'filled','r');
              w =  scatter3(0-dx,1-dy,0-dz,100,'filled','r');
              w =  scatter3(1-dx,0-dy,0-dz,100,'filled','r');
              w =  scatter3(0-dx,0-dy,0-dz,100,'filled','r');
              w =  scatter3(0+dx,0+dy,0+dz,100,'filled','r');
               text(0.2,-0.3,-0.3,'$\omega_4^1$','fontsize',16,'interpreter','latex')
              
        end   
        if i ==3
             % plot3([0 0],[0 0],[1 1+1/15],'r','linewidth',2) ;
             % plot3([1 1],[0 0],[1 1+1/15],'r','linewidth',2) ;
             % plot3([0 0],[1 1],[1 1+1/15],'r','linewidth',2) ;
             % plot3([0 0],[0 0],[0 1/15],'r','linewidth',2) ;  
             % plot3([0+0.4 0+0.4],[0 0],[0.3 0.3+1/15],'r','linewidth',2) ; 
             %  plot3([0+0.4 0+0.4],[0 0],[1+0.3 1+0.3+1/15],'r','linewidth',2) ;
             % plot3([1 1],[0 0],[1 1+1/15],'r','linewidth',2) ;
             % plot3([0+0.4 0+0.4],[1 1],[1+0.3 1+0.3+1/15],'r','linewidth',2) ;
             % plot3([0+0.4 0+0.4],[0 0],[1+0.3 1+0.3+1/15],'r','linewidth',2) ;  
             % plot3([0+0.4 0+0.4],[0 0],[0.3 0.3+1/15],'r','linewidth',2) ; 
              w =  scatter3(0-dx,0-dy,1-dz,40,'filled','r');
              w =  scatter3(0-dx,1-dy,0-dz,40,'filled','r');
              w =  scatter3(1-dx,0-dy,0-dz,40,'filled','r');
              w =  scatter3(0-dx,0-dy,0-dz,60,'filled','r');
              w =  scatter3(0+dx,0+dy,0+dz,40,'filled','r');
              w =  scatter3(1-dx,1-dy,0-dz,30,'filled','r');
              w =  scatter3(0-dx,1-dy,1-dz,30,'filled','r');
              w =  scatter3(1-dx,0-dy,1-dz,30,'filled','r');
              
              w =  scatter3(0+dx,0+dy,1+dz,30,'filled','r');
              w =  scatter3(0+dx,1+dy,0+dz,30,'filled','r');
              w =  scatter3(1+dx,0+dy,0+dz,30,'filled','r');
              w =  scatter3(0+dx,0+dy,0+dz,30,'filled','r');
              w =  scatter3(0+dx,0+dy,0+dz,30,'filled','r');
              text(0.2,-0.3,-0.3,'$\omega_4^2$','fontsize',16,'interpreter','latex')
             
        end   
        if i ==5
            for j = 1:16
                 %plot3([X(j) X(j)],[Y(j) Y(j)],[Z(j) Z(j)+1/16],'r','linewidth',2);
                 w =  scatter3(X,Y,Z,30,'filled','r');
                
            end
          text(0.2,-0.3,-0.3,'$\bar{\omega}_4$','fontsize',16,'interpreter','latex')
        end  
        set(gca,'visible','off')
        view(-20,47)
    end   
  end

  
  annotation('arrow',[0.07 0.07],[0.84 0.2],'linewidth',2)
  annotation('arrow',[0.16 0.8],[0.92 0.92],'linewidth',2)
  
  annotation1 = annotation(...
  'textbox',...
  'Position',[0.4179 0.9095 0.2607 0.07857],...
  'LineStyle','none',...
  'FitHeightToText','off',...
  'FontName','Arial',...
  'FontSize',16,...
  'String',{'$k$ (iteration)'},...
  'Interpreter','latex');
  
 axes('position',[0 0 0.1 1])  
text(0.48,0.3,'$n$ (system size)','fontsize',16,'interpreter','latex','rotation',90)
  set(gca,'visible','off')
 
