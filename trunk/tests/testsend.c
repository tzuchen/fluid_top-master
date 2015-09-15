static char help[] = "A Test for Socket programming.";

#include "petscksp.h"
#include "petscda.h"


int main(int argc,char **args)
{
 PetscErrorCode ierr;
 PetscViewer    socketviewer;
 PetscInt       rank,size,Istart,Iend;
 Vec		x;
 PetscScalar    q = 1403;
 int 		portnumber = 5050,i,s = 10;
 int            fd; 
//PetscViewer_Socket *vmatlab;

  PetscInitialize(&argc,&args,(char *)0,help);
 
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  // Rank 0 connects to socket 5050
     PetscViewerSocketOpen(PETSC_COMM_WORLD,0,portnumber,&socketviewer);    
     PetscViewerBinaryGetDescriptor(socketviewer,&fd); 
    
          // Create a vector x 
 	  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
      	  ierr = VecSetSizes(x,PETSC_DECIDE,s);CHKERRQ(ierr);
      	  ierr = VecSetType(x, VECMPI);CHKERRQ(ierr);
          ierr = VecGetOwnershipRange(x,&Istart,&Iend);CHKERRQ(ierr);
 
          for (i =Istart; i<Iend; i++) { ierr = VecSetValue(x,i,i+30,INSERT_VALUES);}
          ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);          
          VecView(x, 0);
       
	
    

/////////////////////////////////////////////////////////////////////////////////////////

           // Send the vector x
           VecView(x, socketviewer);
          
           // Send a scalar q
           if (rank == 0){                  
	   fprintf(stdout,"q = %f\n",q); 
           PetscBinaryWrite(fd,&q,1,PETSC_DOUBLE,PETSC_FALSE);           
           }

           // Send a Integer s
           if (rank == 0){                  
	   fprintf(stdout,"s = %i\n",s); 
           PetscBinaryWrite(fd,&s,1,PETSC_INT,PETSC_FALSE);           
           }


   ierr = VecDestroy(x);CHKERRQ(ierr);
   ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
