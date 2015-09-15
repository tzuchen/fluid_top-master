% Given A50 and A100
%
% Tzu-Chen 6/1/2006

clear varx h hk
A   = A100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k  = 1e-3;
dt = 1e-3;
p  = 2;

n     =sqrt(size(A,1));
[M]   = rdwalkM(n,p,k,dt);

klist = [35];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dx = 1/n;
    xs = dx/2 :dx:1-dx/2;
    ys = xs;
    [Xs,Ys] = meshgrid(xs,ys);

 
    Xa = 1*sin(2*pi*Ys);
 xa0{1} = reshape(Xa,(n)*(n),1); 

     Xa = 1*sin(2*pi*Xs);
 xa0{2} = reshape(Xa,(n)*(n),1); 
 
    Xa(1:fix(n/2),:) = 1;
    Xa(fix(n/2)+1:end,:) = 0;
 xa0{3} = reshape(Xa,(n)*(n),1);

    Xa(:,1:fix(n/2)) = 0;
    Xa(:,fix(n/2)+1:end) = 1;
 xa0{4} = reshape(Xa,(n)*(n),1);

  nofcase = size(xa0,2);  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:nofcase
  [x,Xa1,varx(i,:)] = testAs(A,M,klist,xa0{i});
 
end

close all
figure
hold on

for i=1:nofcase
  h(i) = semilogy(varx(i,:));
  
end
set(h,'linewidth',2);
set(gca,'YScale','log')
grid on
xlabel('iterations')
ylabel('\sigma')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
