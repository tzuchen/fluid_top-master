
    n = 500
    iter=50;
    mapfunction = @standardmap;
    param{1}= 0.2;

     [Astruc] = maprefine2(n,[],mapfunction,param{1});
     [M]      = fourierdiffuseM(n,1e-8,1,1);  

    dy    = 1/n;
    y     = dy/2:dy:1-dy/2;
    z     = dy/2:dy:1-dy/2;
    [Yi,Zi] = meshgrid(y,z);

 
    xM = 1;
    xm = -1;
    clear X
    
    x  =[];
    xd = [];



for i = 1:iter

    
    %Xa1 = 1*cos(2*pi*Yi);
    %xa1 = reshape(Xa1,n^2,1);
    if i ==1
      Xa1 = 1*cos(2*pi*Yi);
      xa1 = reshape(Xa1,n^2,1);
    else
      xa1  = Astruc.A*xa1;
      Xa1  = reshape(xa1,n,n);
    end
   
    x = [ x xa1]; 
    X{i}  = Xa1;
    [Yi,Zi]=feval(mapfunction,Yi,Zi,1,-1,param);
    
    
        %  imagesc(Xa1,[xm xM]);
        %  axis equal   
        %  axis([1 n 1 n]);
        %  box on
        %  set(gca,'ytick',[])
        %  set(gca,'xtick',[])
        %  set(gca,'ydir','normal');
        %  Mm(i) = getframe;
         
end




for i =1:iter
             Xf     = fft2(X{i}); 
             Xd{i}   = real(ifft2(Xf.*M.M));
             xd      = [xd reshape(Xd{i},n^2,1)];
end


A = abs(x\xd);
A = A';
%A

for i = 1:iter
 A(i,:) = A(i,:)/sum(A(i,:));

end

Hf = diag(ones(iter-1,1),1);
Hf(end,end) = 1;

A = Hf*A;



