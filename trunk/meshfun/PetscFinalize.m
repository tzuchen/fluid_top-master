function err = PetscFinalize(p)
 
 % send a scalar 0, means end!
 send(p,0);
 pause(1)
 if exist('p','var')
  closeport(p); 
 end
 
 err = 0;
