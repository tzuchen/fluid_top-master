function [A,b] = generateAb(mu,L,D,lr,dr,fr,cdinfo)   

t0 = clock;
np = prod(cdinfo.grid);
O  = sparse(np,np);

%A  = -[ mu*L -D' ; 
%       -D     O  ];
%            
%b  = -[ (-mu*lr - fr);
%             dr      ];

% Tzuchen 11-4-2007 make the change, to make u,and p in the same order
A  = -[ L    -D' ; 
       -D     O  ];
            
b  = -[ (-lr - fr/mu);
             dr      ];


if cdinfo.showtime==1
  disp(sprintf('Matlab: Time to generate A and b : %f sec ',etime(clock,t0))); 
end
