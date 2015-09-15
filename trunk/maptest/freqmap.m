function   [P,x0,statelist,X0,Ma] = freqmap(n,cutratio,iter,ic,mapfunction,param)



 k = 0;


 [Astruc] = maprefine2(n,[],mapfunction,param{:}); 

 [Mdiff]  = fourierdiffuseM(n,k,1,1);

 nlist =  [n/2:-1:1,0:1:n/2-1];

 P     =  ones(n,1)*nlist.^2+ (ones(n,1)*nlist.^2)';
 M     =  zeros(n,n);
 M0    =  M;

 dx = 1/n;
 xs = dx/2 :dx:1-dx/2;
 ys = xs;
 [Xs,Ys] = meshgrid(xs,ys);
 Xa  = zeros(n,n);

 % ic
  if ic == 'cosx'
      Xa = 1*cos(2*pi*Ys);
  elseif ic == 'cosy'
      Xa = 1*cos(2*pi*Xs);
  elseif ic == 'squx'
      Xa(1:fix(n/2),:) = 0.5;
      Xa(fix(n/2)+1:end,:) = -0.5;
  elseif ic == 'squy'
      Xa(:,1:fix(n/2)) = 0.5;
      Xa(:,fix(n/2)+1:end) = -0.5;
  elseif ic == 'rand'
      Xa = rand(n,n);
  end

 M  = fft2(Xa);

 M  = M.*abs(M>1e-5); 
 M  = M/max(max(abs(M)));
 statelist    = [1 ;find(M)];
 icf          =  reshape(M,n^2,1);
 normf        =  [];

X     = ifft2(M);
x     = reshape(X,n^2,1);
 for i = 1:iter
    %X     = ifft2(M);
    %x     = reshape(X,n^2,1);
    if i>1
       x     = Astruc.A*x; 
       X     = reshape(x,n,n);
    else
       X0    = X;
    end 
    normf = [normf sqrt(var(x))];  

    M     = fft2(X);
   Ma  = M;
 %Mp = M0;
 %Mp(statelist) = 1;
 %spy(fftshift(Mp));
 %Mn(i) = getframe;
 


    %[iadd,jadd] = find(abs(M)>cutratio);
    stateadd =  find(abs(M)>cutratio);
    %stateadd =stateadd(find(jadd>=n/2-1));
    statelist = union(statelist,stateadd);
        
 end

    ns  = length(statelist);
    P   = sparse(2*ns,2*ns);
    P11 = zeros(ns,ns);
    P12 = zeros(ns,ns);
    length(statelist); 

fftw('planner','patient');
tic
for i = 1:ns

    M = M0;
    M(statelist(i))= 1;

    X     = ifft2(M);
    x     = reshape(X,n^2,1); 

    x     = Astruc.A*x; 
    X     = reshape(x,n,n);
    M     = fft2(X); 

    %M     = M.*Mdiff.M; 
    M     = M.*(abs(M)>=cutratio);

    Ms    = M(statelist);

    P11(i,:)   =  (reshape(real(Ms),1,ns));
    P12(i,:)   =  (reshape(imag(Ms),1,ns));
 
end
toc


%tic
% Fp = my2dDFT(n,statelist,1);
%toc
%%tic
%% Fn = my2dDFT(n,statelist,-1);
%%toc
%tic
% Ap = (Fp*Astruc.A*Fp')./n^2;
%toc


  P = sparse([P11 P12;-P12 P11]);


  normf = normf/normf(1);  


normx = [];
sigx  = [];


x= icf(statelist);
y = x;
x = [real(x);imag(x)];
x0 = x;



for i= 1:iter
 
   sigx = [sigx norm(x(2:end))];
   x = P'*x;

end

close all
figure


sigx = sigx/sigx(1);
subplot(2,1,1)
hold on
h= plot(normf,'r');
set(h,'linewidth',2);
h= plot(sigx,'g');
set(h,'linewidth',2);
grid on
box on
hold off

%figure

subplot(2,1,2)
hold on
h = semilogy(normf,'r');
set(h,'linewidth',2);
h = semilogy(sigx,'g');
set(h,'linewidth',2);
set(gca,'YScale','log')
grid on
box on
hold off














