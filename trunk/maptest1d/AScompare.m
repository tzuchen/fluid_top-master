function AScompare(n,mun,S)

%  Compare the numerical A
%  and the analytical P
%
%  Tzu-Chen Liang 6-7-2007
close all
figure 
hold on



dx = 1/n;
x = dx/2:dx:1-dx/2;
[A,B] = maprefine1d3(n,S);
y0 = zeros(n,1);
y0(1:mun*n) = 1/mun;
y = y0;

for k = 1:3

   w = wfunS(x,k,mun,S) 
   y = A'*y;
 
   h1 = plot(x,w)
   h2 = plot(x,y,'r--')

   

end  

axis([0 1 0 10])
