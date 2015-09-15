T= 1;
N= 500;
dt =T/N;
close all
figure 
hold on
clear W
clear ito

for i = 1:1000

dW = sqrt(dt)*randn(1,N);
W(i,:)  = cumsum(dW);

%plot(W(i,:))

ito2(i,:) = cumsum([0,W(i,1:end-1)].*dW);
ito3(i,:) = cumsum(([0,W(i,1:end-1)].^2).*dW);

end

tl = 0:dt:1;



plot((tl.^2)/2)
plot(var(ito2))

plot((tl.^3),'r')
plot(var(ito3),'r')
