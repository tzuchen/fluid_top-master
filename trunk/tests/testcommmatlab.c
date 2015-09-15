static char help[] = "A Test for Socket programming.";

#include "petscksp.h"
#include "petscda.h"
#include <stdlib.h>

extern const int VecReceive(PetscViewer, VecType,Vec*);
extern const int MatReceive(PetscViewer, MatType,Mat*);


int main(int argc,char **args)
{
 PetscErrorCode ierr;
 PetscViewer    socketviewer;
 PetscInt       rank,size;
 Vec		x,y;
 Mat            A;
 PetscScalar    *a;
 int 		vsize,portnumber = 5050;
 
//PetscViewer_Socket *vmatlab;

  PetscInitialize(&argc,&args,(char *)0,help);
 
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
  // Rank 0 connects to socket (use PETSC_DEFAULT, not necessary to be #5050)
     //PetscViewerSocketOpen(PETSC_COMM_WORLD,0,portnumber,&socketviewer); 
     PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer); 

  // Example for receiving and sending a Scalar
     // It becomes a VECSEQ with size one.
     VecReceive(socketviewer,PETSC_NULL,&x);
     //VecView(x,0);     
     // We cannot send VECSEQ type vector, so use PetscRealView, and notice the a is a pointer.
     VecGetArray(x,&a);
     PetscRealView(1,a, socketviewer);
     

  // Example for receving and sending a SEQ Vector
     // If the size of the receiving vector is larger than 1, it becomes a VECMPI type vector. 
    
     VecReceive(socketviewer,VECSEQ,&x);
     VecGetArray(x,&a);
     VecGetSize(x,&vsize); 
     PetscRealView(vsize,a, socketviewer);

  // Example for receving and sending a MPI Vector
  
     VecReceive(socketviewer,VECMPI,&y);
     //VecView(y,0);
     // Use VecView to send a VECMPI type vector 
     VecView(y, socketviewer);

  // Example for receving and sending a MPI Matrix
    
      MatReceiveTranspose(socketviewer,MATMPIAIJ,&A);
     //MatView(A,0);
     // MatView(A,socketviewer);



/////////////////////////////////////////////////////////////////////////////////////////
   ierr = VecDestroy(x);CHKERRQ(ierr);
   ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

