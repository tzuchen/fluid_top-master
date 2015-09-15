function []=mapplot(n,r1,r2,K)

%close all 
hold on

h0(1)=plot([0,0],[-2*n,2*n]);
h0(2)= plot([-2*n,2*n],[0,0]);

h1 = plot(r1,r2,'ro');
set(h1,'markersize',5)


h2(1) = plot([n,n],[n,-n]);
h2(2) = plot([n,-n],[-n,-n]);
h2(3) = plot([-n,-n],[-n,n]);
h2(4) = plot([-n,n],[n,n]);
set(h2,'linestyle','--','linewidth',2);


%h3 = plot([-2*n,2*n],[r1+r2,r1+r2]);
%set(h3,'color','r')


w1 = [-K*(r1+r2)+r1,r1+r2];
w2 = [ K*(r1+r2)+r1,r1+r2];


% make w1<w2
[w1,w2] = myswitch(w1,w2);


h4(1) = plot(w1(1),w1(2),'xr');
h4(2) = plot(w2(1),w2(2),'xr');

if and(r1+r2<=n,r1+r2>=-n)
  w3 = [-n,r1+r2];
  w4 = [n ,r1+r2];
  h4(3) = plot(w3(1),w3(2),'xr');
  h4(4) = plot(w4(1),w4(2),'xr');
else 
  w3 = [];
  w4 = [];
end





h5 = plot([w1(1),w2(1)],[w1(2),w2(2)]);
set(h5,'color','g')


pw1(1) = -K*(w1(1)+w1(2))+w1(1);
pw1(2) =  K*(w1(1)+w1(2))+w1(1);
[pw1(1),pw1(2)] = myswitch(pw1(1),pw1(2));

h6(1) = plot([pw1(1),pw1(2)],[w1(1)+w1(2),w1(1)+w1(2)]);

pw2(1) = -K*(w2(1)+w2(2))+w2(1);
pw2(2) =  K*(w2(1)+w2(2))+w2(1);
[pw2(1),pw2(2)] = myswitch(pw2(1),pw2(2));
h6(2) = plot([pw2(1),pw2(2)],[w2(1)+w2(2),w2(1)+w2(2)]);

if isempty(w3)~=1
  pw3(1) = -K*(w3(1)+w3(2))+w3(1);
  pw3(2) =  K*(w3(1)+w3(2))+w3(1);
  [pw3(1),pw3(2)] = myswitch(pw3(1),pw3(2));
  h6(3) = plot([pw3(1),pw3(2)],[w3(1)+w3(2),w3(1)+w3(2)]);

  if w3(1)>w1(1)
    fill([pw1(1) pw1(2) pw3(2) pw3(1)],[w1(1)+w1(2),w1(1)+w1(2),w3(1)+w3(2),w3(1)+w3(2)],'r' )
  end
end

if isempty(w4)~=1
  pw4(1) = -K*(w4(1)+w4(2))+w4(1);
  pw4(2) =  K*(w4(1)+w4(2))+w4(1);
  [pw4(1),pw4(2)] = myswitch(pw4(1),pw4(2));
  h6(4) = plot([pw4(1),pw4(2)],[w4(1)+w4(2),w4(1)+w4(2)]);

  if w2(1)>w4(1)
    fill([pw2(1) pw2(2) pw4(2) pw4(1)],[w2(1)+w2(2),w2(1)+w2(2),w4(1)+w4(2),w4(1)+w4(2)],'r' )
  end
end


if isempty(w4)==1
   fill([pw1(1) pw1(2) pw2(2) pw2(1)],[w1(1)+w1(2),w1(1)+w1(2),w2(1)+w2(2),w2(1)+w2(2)],'r' )
end





axis equal
%axis([-2*n 2*n -2*n 2*n])
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%
function [x1,x2] = myswitch(x1,x2);

  if x1(1)>x2(1)
     xt = x2;
     x2 = x1;
     x1 = xt;

  end










