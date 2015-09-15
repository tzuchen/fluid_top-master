function err = PetscFinalize(p)
 
 % send a scalar, means end!
 send(p,1);

 if exist('p','var')
 closeport(p); 
 end
 endflag = 1;
 save('endflag.mat','endflag');  

 if exist('taskflag.mat','file')
   delete taskflag.mat
 end
 err = 0;
