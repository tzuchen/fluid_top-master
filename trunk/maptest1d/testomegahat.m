function testomegahat()

close all
klist = 0:9;
nlist = 1:6;
dx  = 1/1000;
y= dx/2:dx:1-dx/2;
colorlist = 'gcmybgcmybgcmyb'
for k = klist
  for n = nlist
   
    if k<n+3
    w = 0.5^k*Tp(y*0.5^k)./T(mun(n));
    if k<n 
      w(find(y>munk(n,k)))=0;
    end

    if k ==n+2
    stairs(y,w,colorlist(n));
    end

    hold on
   axis([0 1 0 1])
   end
 end
end

figure
hold on
plot(y,Tp(y)./T(1/2))
plot(y,0.5*Tp(0.5*y)./T(1/4),'r')
plot(y,T(y))

%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = T(y)
x = sin(pi*y/2).^2;

function x = Tp(y)
x = 0.5*pi*sin(pi*y);

function mu = mun(n)
 mu = 2^(-n+1);

function l = munk(n,k)
 l = 2^(k-n+1);


