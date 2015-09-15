function [v,lan,Dv,Dlan] = mypower(A,plist)
% This is the verson that works for large sparse problem.
 n           = size(A,1);
 v           = zeros(n,1);
 vk          = rand(n,1);
 vkn         = vk / norm(vk);
 
 p           = length(plist);
 Dvkn        = sparse(n,p);
 DA          = sparse(n,p);
 M           = sparse(n,p);
 itable      = mod(plist,n);
 itable(find(itable==0)) = n;
 jtable      = fix((plist-1)./n);
 jtable      = jtable + 1;

 
  
 
 while(norm(vkn-v)>1e-10)
   v         = vkn;   
   Dv        = Dvkn;
   
   vk        = A*v;
   for i =1:p 
        M(itable(i),i) = v(jtable(i));        
   end
   Dvk       = M + A*Dv;  
    
   
   vknorm    = norm(vk); 
   vkn       = vk / vknorm;
   Dvknorm   = (Dvk'*vk)/vknorm;
   Dvkn      = Dvk/vknorm-vk*Dvknorm'/vknorm^2;  
 end
 
 
 v = vkn;
 Dv = Dvkn;  
  
 lan = (A*v)./v; 
 lan = lan(1);
 
 a1 = A(1,:);
 
 v1  = v(1);
 Dv1 = Dv(1,:);
 
 vexp = zeros(n,1);
 vexp(1) = 1;
 
 
     for i =1:p 
        DA(jtable(i),i) = vexp(itable(i));        
     end
     
 Dlan = (v'*DA+ a1*Dv- lan*Dv1)/v1; 
 
 
 
