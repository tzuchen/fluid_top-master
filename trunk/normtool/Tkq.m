function t = Tkq(k,q,eps)


 %M =[1 1;0 1];
 M = [2 1;1 1];
 
 %K = 0.5*2*pi;
 K = 1e-3;

 Q = (k'*M-q')';

 if Q(2)==0
 del = 1;
 else
 del = 0;
 end

 %t = exp(-eps*(q'*q))
 t = del*(i.^Q(1,:)).*besselj(Q(1,:),K*(k(1,:)+k(2,:))); 
