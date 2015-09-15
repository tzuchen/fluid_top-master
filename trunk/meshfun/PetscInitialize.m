function err = PetscInitialize(np,filetoexecute,option)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Set options
% default value 
  rtol     = 1e-4;
  ksp_type = 'symmlq';

  
 %if nargin == 3
 %   if isfield(option,'rtol')
 %       rtol   = option.rtol;
 %   end
 %   if isfield(option,'ksp_type')
 %       ksp_type   = option.ksp_type;
 %   end
 %
 %else 
 %   rtol = 1e-5;
 %   ksp_type = 'symmlq'; 
 %end
 
if nargin == 3 
 F = fieldnames(option);
 n = size(F,1); 
 for i = 1:size(F,1)
   opt{i} = [' -',F{i},' ',num2str(getfield(option,F{i})) ];
 end
else
  n = 0;
end  
  opt{n+1} = [' -rtol ', num2str(rtol) ];
  opt{n+2} = [' -ksp_type ', ksp_type ];
  opt{n+3} = [' --hostfile mpd.hosts'];
  opt{n+4} = [' -withMatlab 1']; 

    %opt{1}  = [' -rtol ', num2str(rtol) ]; 
    %opt{2}  = [' -ksp_type ', ksp_type ];    
    %opt{3}  = [' -log_summary '];
    optset  = [opt{:}];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%err = launch(['./fluidsolver' optset],np);
 launch(['./',filetoexecute,optset],np);
err=0;

if err==0
   disp(sprintf('Petsc solver launched successfully!'))
else
   disp(sprintf('Failed to launch Petsc solver! Please stop the program.'))
end



