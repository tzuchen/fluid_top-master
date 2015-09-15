% This demo shows how eigenvalue changes when grid size gets smaller  
%

%testmeshfun
%close all
%figure 
hold on
clear A dA mu D

for n = 20:5:55

  [A{n},dA{n},mu(n),D{n}] = testmarkovmap(n,n,cdinfo,sol);
  h = plot(sort(diag(abs(D{n}))),'r');
  set(h,'linewidth',2);


end

grid on
xlabel('Number')
ylabel('Eigenvalue')

  
