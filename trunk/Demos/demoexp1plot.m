close all
figure 
hold on

X = [0.3 0.7 0.7 0.3 0.3]
Y = [0.3 0.3 0.7 0.7 0.3]
fill(X,Y,'g')

h(1) = plot([0.3 0.3],[0.3,0.7])
h(2) = plot([0.7 0.7],[0.3,0.7])
h(3) = plot([0.3 0.7],[0.3,0.3])
h(4) = plot([0.3 0.7],[0.7,0.7])

set(h,'linewidth',2)
axis equal
axis([0,1,0,1])

box on

