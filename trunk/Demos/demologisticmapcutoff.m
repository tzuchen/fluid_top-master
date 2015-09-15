% Demo logisticmap cutoff
%
%
%
% Tzu-Chen Liang  4-28-2007
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nlist = [1e2 1e3 1e4 1e5 1e6 ]
dl = length(nlist);
nit = 30;
m1all = [];
m2all = [];
m3all = [];
m1all2 = [];
m2all2 = [];
m3all2 = [];


for i= 1:dl
A = maprefine1d3(nlist(i),@logisticmap);
[X,m1,m2,m3,x0]=Atest1d(A,'cosx',nit,[],0);
m1all = [m1all;m1 ];
m2all = [m2all;m2 ];
m3all = [m3all;m3 ];

% First find the invariant distribution by simulating 10*nit
%[X,m1,m2,m3,x0]=Atest1d(A','sing',10*nit);
[X,m1,m2,m3,x0]=Atest1d(A','sing',nit,[],0);
m1all2 = [m1all2;m1 ];
m2all2 = [m2all2;m2 ];
m3all2 = [m3all2;m3 ];

end

for i = 1:dl
  w{i} = sprintf('n = %d',nlist(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
hold on

h = plot(1:nit,m3all);
set(h,'linewidth',2);
legend(w)
grid on
box on
xlabel('iteration')
ylabel('var(f^k)')
title('cosine function advected by A')
axis([1 nit 0 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure 
hold on

h = plot(1:nit,m2all2);
set(h,'linewidth',2);
legend(w)
grid on
box on
xlabel('iteration')
ylabel('|x-\pi|/2')
title('x=[1 0 0 ...]^T evolved by A^T')
axis([1 nit 0 1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




