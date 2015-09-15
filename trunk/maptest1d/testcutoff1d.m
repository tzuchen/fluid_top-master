
%nlist = [1e2 1e3 1e4 1e5 1e6 1e7];
nlist = [1e2 1e3 1e4 1e5 1e6 ];
m1l=[];
m2l=[];
m3l=[];



for i = 1:length(nlist)
   n = nlist(i);
 
  A = maprefine1d3(n,@logisticmap);
  %[X,m1,m2,m3,x0]=Atest1d(A','sing',1000);
  %[X,m1,m2,m3,x0]=Atest1d(A','sing',100,X);
  %[X,m1,m2,m3,x0]=Atest1d(A','sing',100);
  [X,m1,m2,m3,x0]=Atest1d(A,'cosx',30,[],0);
  %[X,m1,m2,m3,x0]=Atest1d(A','sing',30,[],0);

  m1l = [m1l; m1];
  m2l = [m2l; m2];
  m3l = [m3l; m3];
end


close all
figure
plot(m1l');
figure
plot(m2l');
figure
plot(m3l');
