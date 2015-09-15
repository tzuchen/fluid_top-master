
n      = 100;
dx     = 1/n;
x      = dx/2:dx:1-dx/2;
y      = dx/2:dx:1-dx/2;
slope  = [];
varmat = [];

[X,Y] = meshgrid(x,y);

freq = [1 2 3 4];

for i = 1:length(freq)
   
    Z           = sin(2*pi*X*freq(i));
    var2        = diffsimu(100,Z,1e-3,20);
    var2        = var2/var2(1);
    varmat(i,:) = var2;
    slope(i)    = (log(var2(end))-log(var2(end-10+1)))/10;
end

semilogy(varmat')
%plot(freq,slope)

slope(2)/slope(1)
slope(3)/slope(1)
slope(4)/slope(1)

-slope(1)/4/pi^2/1e-3
