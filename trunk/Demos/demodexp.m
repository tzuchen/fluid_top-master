% Demos for two kinds of cutoff
%
% 
% Tzu-Chen Liang

close all
clear f

  k = 1:0.1:10;

close all
figure 
hold on

for n = 1:10
   f(n,:)  =exp(-exp(-n+k-3));
end

h = plot(k,f','b');
set(h,'linewidth',2)
grid on
box on
xlabel('k')
ylabel('D')


figure 
hold on
k = -100:0.1:100
for n = 1:10-2
   k1 = k/(log(-log(0.5))+3+n);  
   h = plot(k1,exp(-exp(-n+k-3)));
   set(h,'linewidth',2,'color','b');
   
end
grid on
box on
axis([0 2 0 1])
xlabel('normalized k')
ylabel('D')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear f

  k = 1:0.1:5;


figure 
hold on

for n = 1:5
   f(n,:)  =exp(-(-n+k-3));
end

h = plot(k,f','b');
set(h,'linewidth',2)
grid on
box on
xlabel('k')
ylabel('D')


figure 
hold on
k = -100:0.1:100
for n = 1:10
   k1 = k/(-log((0.5))+3+n)  ;
   h = plot(k1,exp(-(-n+k-3)));
   set(h,'linewidth',2,'color','b');
   
end
grid on
box on
axis([0 2 0 1])
xlabel('normalized k')
ylabel('D')





