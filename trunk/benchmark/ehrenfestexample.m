
iter = 15000;
dlist = [ 100 200 400 1000 2000];
%dlist= [10000]
dl    = length(dlist);

close all
figure 
hold on
nlist = []

for i=1:dl
 
  n = dlist(i)
  [A,xe] = ehrenfestA(n);
  x0 = zeros(n+1,1);
  xe = binopdf(0:n,n,0.5)';
  x0(1) = 1;
  x = x0;
  nlisti = [];
  for j = 1:iter
     x = A'*x;
     nlisti = [nlisti sum(abs(x-xe))/2];
     %nlist = [nlisti sqrt(sum(xe.*(x./xe).^2))];
  end  
   nlist(i,:) =nlisti;

  


end

h= plot(nlist')
set(h,'linewidth',2)

for i = 1:dl
  w{i} = sprintf('n = %d',dlist(i));
end
legend(w);



grid on
box on
xlabel('Iterations')
ylabel('total variation distance')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


noval = 0.5

for i = 1:dl
  for j = 1:iter-1
    if nlist(i,j+1)>=nlist(i,j);
      nlist(i,j+1) = nlist(i,j)-1e-9; 
    end
  end
end


for i = 1:dl
  itern(i)   = interp1(nlist(i,:),0:iter-1,noval);
end

clist = 'bgrcmyk'

figure
hold on

for i = 1:dl

h(i) = plot([0:iter-1]/itern(i),nlist(i,:));
set(h(i),'color',clist(i));  
end
legend(w);
axis([0 5 0 1])
set(h,'linewidth',2);

box on
grid on
xlabel('nornalized iteration')
ylabel('total variation distance')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear xi yi areai
figure 
hold on
box on
xi = 0:1e-4:3;
for i = 1:dl

 yi(i,:) = interp1([0:iter-1]/itern(i),nlist(i,:),xi,'spline');
 areai(i) = trapz(xi,abs(yi(i,:)-0.5));
end

g1 = plot(log(dlist.*log(dlist)/4),3-areai)
set(g1,'linewidth',2);
g2 = plot(log(dlist.*log(dlist)/4),3-areai,'s')
set(g2,'markersize',10,'linewidth',2);





