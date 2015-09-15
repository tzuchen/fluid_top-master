% Demo for the polyintersect function

close all
clear X1 X2
figure
hold on
color = 'rgby';
count = 0;

%ny = 20;
%nz = ny;
%dy = 1/ny;
%[Y0,Z0]= meshgrid([0:dy:1],[0:dy:1]);
%[Y,Z]  = standardmap(Y0,Z0,1,1);
load y0filew
load z0filew
ny = 51; nz = 51;
dy = 1/(ny-1);
Y0 = reshape(y0filew,ny,nz);
Z0 = reshape(z0filew,ny,nz);
Y =  reshape(yefilew,ny,nz);
Z =  reshape(zefilew,ny,nz);





A = sparse((ny-1)*(nz-1),(ny-1)*(nz-1));
debugmode = 0;

for  i = 1:ny-1
  for j = 1:nz-1;
     for k = 1:ny-1;
         for l = 1:nz -1

  	   if abs(i-k)<=8
    	      if abs(j-l)<=8 

                
		X1(1,:) = reshape(Y0(i:i+1,j:j+1),1,4);
	        X1(2,:) = reshape(Z0(i:i+1,j:j+1),1,4);
		X2(1,:) = reshape( Y(k:k+1,l:l+1),1,4);
	        X2(2,:) = reshape( Z(k:k+1,l:l+1),1,4);

		X1(:,:) = X1(:,[1 2 4 3]);
		X2(:,:) = X2(:,[1 2 4 3]);
            
                x1xm = min(X1(1,:));
                x1xM = max(X1(1,:));
                x2xm = min(X2(1,:));
                x2xM = max(X2(1,:));
                x1ym = min(X1(2,:));
                x1yM = max(X1(2,:));
                x2ym = min(X2(2,:));
                x2yM = max(X2(2,:));



                if and(~or(and(x1xm<x2xm,x1xM<x2xm),...
                           and(x2xm<x1xm,x2xM<x1xm)),...
                       ~or(and(x1ym<x2ym,x1yM<x2ym),...
                           and(x2ym<x1ym,x2yM<x1ym)) )
                                    

                   count = count+1;
		   [X3, area] = polyintersect(X1,X2);
                   A((nz-1)*(j-1)+i,(nz-1)*(l-1)+k) = area/dy/dy;
 
              if debugmode ==1
                    close all
                    figure
                    hold on
                    patch(X1(1,:),X1(2,:),'r')
                    patch(X2(1,:),X2(2,:),'g')
                    patch(X3(1,:),X3(2,:),'b')
                    axis equal
                    pause
             else    
                    patch(X3(1,:),X3(2,:),color(mod(count,4)+1))
             end
                   disp(sprintf('%d a intersection area =  %d',count,area))

                end 
	    end
	  end
	end
     end 
  end
end


X1 = zeros((ny-1),(nz-1));
X1((ny-1)/2+1:end,:) = 1;
x1 = reshape(X1,(ny-1)*(nz-1),1);

%figure
%contourf(reshape((A^1)*x1,ny-1,nz-1))
%figure
%contourf(reshape((A^10)*x1,ny-1,nz-1))
%figure
%contourf(reshape((A^100)*x1,ny-1,nz-1)) 


%patch(X1(1,:),X1(2,:),'g')
%patch(X2(1,:),X2(2,:),'r')





