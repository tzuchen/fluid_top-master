close all
 figure
 hold on
 clear w

  eps  = 0.1;
  r    = 0.9;
  d(1) = eps+0.5;



for n = 12:12
   %n    = 10;
   for i = 2:n
     d(i)= (d(i-1)-0.5)*r+0.5;
   end
   dn = 1-d;

   for i = 0:2^n-1
      s = dec2bin(i,n);
      t =  str2num(s')';
      a = prod(d(find(t)));
      b= prod(dn(find(t-1)));
      w(i+1) = a*b ;

    end

stairs([1:2^n]/2^n ,sort(w)*2^n);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% binomial
p = ceil(1/(1-r))
y = binopdf(0:p,p,d(1));

wp =[];
for i = 0:p
     pf = nchoosek(p,i);
   wp = [wp y(i+1)/pf*ones(1,pf) ];

end
size(wp)
y
stairs([1:2^p]/2^p ,sort(wp)*2^p,'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n = 12
k = 9;
y = zeros(1,n+1)';
% for w = 0:n
%     for i=0:w
%         for s=0:n-w
%           
y(w+1)=y(w+1)+((1-2*(i+s)/(n+1))^k)*((-1)^i)*(nchoosek(w,i))*(nchoosek(n-w,s));
%         end
%     end
% end
%     y = y/(2^n);
%wp =[];
y(1)=  1;
A = ehrenfestA(n)
for i =1:k
    y = A'*y;
end


wp =[];
for i = 0:n
     pf = nchoosek(n,i);
   wp = [wp y(i+1)/pf*ones(1,pf) ];

end
sum(wp)
sum(y)
stairs([1:2^n]/2^n ,sort(wp)*2^p,'g');
