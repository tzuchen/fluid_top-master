/* Program usage:  mpirun -np <procs> fluidsolver [-help] [all PETSc options] */ 

static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
  -random_exact_sol		   	 : use a random exact solution vector\n\
  -view_exact_sol   		   	 : write exact solution vector to stdout\n\
  -rtol <error tolerance>          	 : error tolerance\n\   ";




#include "petscksp.h"
#include "petscda.h"
#include <sys/stat.h>
#include <time.h>

extern const int VecReceive(PetscViewer, VecType, Vec*);
extern const int MatReceiveTranspose(PetscViewer, MatType, Mat*);
extern const int ScalarReceive(PetscViewer, PetscScalar**);
extern const int IntReceive(PetscViewer, PetscInt**);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            b,x;   			
  Mat            A,I;        		
  KSP            ksp;    	
  PetscInt       rank,size,*EndTag,its;
  PetscErrorCode ierr;
  PetscScalar    rtol=1.e-12;  
  PetscViewer    socketviewer;
  PetscLogDouble t1,t2,elapsed_time;
  
  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
  // Rank 0 connects to socket, use default socket
     PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer); 
    //PetscViewerSocketOpen(PETSC_COMM_WORLD,0,5050,&socketviewer); 
     ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: socket opened! \n");CHKERRQ(ierr);     
  
  // List of parameters
  //   ierr = PetscOptionsGetScalar(PETSC_NULL,"-rtol",&rtol,PETSC_NULL);CHKERRQ(ierr);

  
  
  // Set ksp
  	  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  	  ierr = KSPSetTolerances(ksp,rtol,1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  	  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);


          ierr = IntReceive(socketviewer,&EndTag);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Problem Number : %d\n",*EndTag);CHKERRQ(ierr);
   while(*EndTag != 0) {


/////////////////////////////// 
/*  Receive A,I and b          */
///////////////////////////////


  // Receive matrix A
    ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&A);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix A received.\n");CHKERRQ(ierr);  
 
  // Receive matrix I
    ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&I);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix I received.\n");CHKERRQ(ierr);  

  // Receive RHS vector b
    ierr = VecReceive(socketviewer,VECMPI,&b);CHKERRQ(ierr);
    //ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: vector b received.\n");CHKERRQ(ierr);   
 
  if (*EndTag==1){
  // Form x by copying b       
    ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
     }      
/////////////////////////////// 
/*  Solver Ax =b             */
///////////////////////////////


  ierr = PetscGetTime(&t1);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,I,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
  ierr = PetscGetTime(&t2);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);CHKERRQ(ierr);
  elapsed_time = t2 - t1;
  ierr = KSPGetIterationNumber(ksp,&its);CHKERRQ(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,"PETSC: x is solved in %4.2e secs with %d iterations \n",elapsed_time,its);
  
  

/////////////////////////////// 
/*  Send x to Matlab socket  */
///////////////////////////////

   
    VecView(x,  socketviewer);
    
    PetscPrintf(PETSC_COMM_WORLD,"PETSC: x has been sent to Matlab socket.\n\n");
     
   
    // Receive EngTag (if EngTag=0, break)
    ierr = IntReceive(socketviewer,&EndTag);CHKERRQ(ierr);
    if (*EndTag < 1) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Received EngTag. Finalize PETSC. \n");CHKERRQ(ierr);
    }else{
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Problem Number : %d\n",*EndTag);CHKERRQ(ierr);
         }

}
/////////////////////////////// 
/*  Destroy Everything       */
///////////////////////////////

  ierr = KSPDestroy(ksp);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}



