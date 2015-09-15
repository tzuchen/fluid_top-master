function mapsimu1d(mapfunction1d)

x = rand;

close 
figure 
hold on

for i = 1:10000
x = feval(mapfunction1d,x);
plot(i,x)
end

