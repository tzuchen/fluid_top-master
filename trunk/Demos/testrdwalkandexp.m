 close all
 figure
 hold on
 clear w

  n = 30;
  k = 200;%1/4*n*log(n);  


  r = exp(-2/n);
  eps = sqrt( n*(1-r)/8 );

  eps = eps*r^k;
 

 % eps  = 0.1;
 % r    = 0.9;
  d(1) = eps+0.5;



for np = 12
   %n    = 10;
   for i = 2:np
     d(i)= (d(i-1)-0.5)*r+0.5;
   end
   dn = 1-d;

   for i = 0:2^np-1
      s = dec2bin(i,np);
      t =  str2num(s')';
      a = prod(d(find(t)));
      b= prod(dn(find(t-1)));
      w(i+1) = a*b ;

    end

stairs([1:2^np]/2^np ,sort(w)*2^np);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%binomial
%p = ceil(1/(1-r))
%y = binopdf(0:p,p,d(1));

%wp =[];
%for i = 0:p
%     pf = nchoosek(p,i);
%     wp = [wp y(i+1)/pf*ones(1,pf) ];
%end
%stairs([1:2^p]/2^p ,sort(wp)*2^p,'r');


%binomial
p = ceil(1/(1-r))
y = binopdf(0:p,p,d(1));


nu = [];
z = [];
for i = 0:p
     pf = nchoosek(p,i);
     z(i+1)= y(i+1)/pf;
     nu(i+1) = pf; 
end
[z,I]= sort(z);

 nu  = [0 cumsum(nu(I)/2^p)];
z = z*2^p;

nu
for i = 1:p
  plot([nu(i),nu(i+1)], [z(i),z(i)],'r');
  plot([nu(i+1),nu(i+1)], [z(i),z(i+1)],'r');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%y = zeros(1,n+1)';
%y(1)=  1;
%A = ehrenfestA(n);
%for i =1:k
%    y = A'*y;
%end

%wp =[];
%for i = 0:n
%     pf = nchoosek(n,i);
 %  wp = [wp y(i+1)/pf*ones(1,pf) ];

%end

%stairs([1:2^(n)]/2^n ,sort(wp)*2^n,'g');


y = zeros(1,n+1)';
y(1)=  1;
A = ehrenfestA(n);
for i =1:k
    y = A'*y;
end

nu = [];
z = [];
for i = 0:n
     pf = nchoosek(n,i);
     z(i+1)= y(i+1)/pf;
     nu(i+1) = pf; 
end
[z,I]= sort(z);

 nu  = [0 cumsum(nu(I)/2^n)];
z = z*2^n;

nu
for i = 1:n
  plot([nu(i),nu(i+1)], [z(i),z(i)],'g');
  plot([nu(i+1),nu(i+1)], [z(i),z(i+1)],'g');
end





