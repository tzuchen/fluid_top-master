
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
#include "MeshStruc.h"
#include "PetscReceive.h"

     
#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            b,c,x,y,adjvec,gradf;   			
  Mat            A,Ac,Ap,Prt,alphamat;        		
  KSP            ksp1,ksp2;    	
  PetscInt       m,n,rank,size, *EndTag,its;
  PetscErrorCode ierr;
  PetscScalar    rtol=1.e-5,fobj;  
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
     ierr = PetscOptionsGetScalar(PETSC_NULL,"-rtol",&rtol,PETSC_NULL);CHKERRQ(ierr);

  
/////////////////////////////// 
/*  Receive A and b          */
///////////////////////////////


  // Receive matrix A
    ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&A);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix A received.\n");CHKERRQ(ierr);  
  // Receive the preconditioner matrix Ac
    ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&Ac);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix Ac received.\n");CHKERRQ(ierr);  
  // Receive the mapping matrix Prt (the transepose of Pr)
    ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&Prt);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix Prt received.\n");CHKERRQ(ierr);  
  // Receive RHS vector b
    ierr = VecReceive(socketviewer,VECMPI,&b);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: vector b received.\n");CHKERRQ(ierr);   
  // Receive RHS vector c 
    ierr = VecReceive(socketviewer,VECMPI,&c);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: vector c received.\n");CHKERRQ(ierr);  

  // Form x,y and adjvec by copying b       
          ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
          ierr = VecDuplicate(b,&y);CHKERRQ(ierr);
          ierr = VecDuplicate(b,&adjvec);CHKERRQ(ierr);

       
  // Form Ap by copying A    
          ierr = MatDuplicate(A,MAT_COPY_VALUES,&Ap);CHKERRQ(ierr);
  // Set ksp
  	  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp1);CHKERRQ(ierr);
  	  ierr = KSPSetTolerances(ksp1,rtol,1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  	  ierr = KSPSetFromOptions(ksp1);CHKERRQ(ierr);

  	  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp2);CHKERRQ(ierr);
          ierr = KSPSetTolerances(ksp2,rtol,1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  	  ierr = KSPSetFromOptions(ksp2);CHKERRQ(ierr);
  

  // Get some number 
          // m: number of x
          // n: number of alpha
          ierr = MatGetSize(Prt,&m,&n);CHKERRQ(ierr);

  // Set vector gradf 
	  ierr = VecCreate(PETSC_COMM_WORLD,&gradf);CHKERRQ(ierr);
      	  ierr = VecSetSizes(gradf,PETSC_DECIDE,m);CHKERRQ(ierr);
      	  ierr = VecSetType(gradf, VECMPI);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Completed data receiving and setting!\n\n");CHKERRQ(ierr);

 /////////////////////////////// 
/*  Receive alphavec         */
/////////////////////////////// 
          ierr = IntReceive(socketviewer,&EndTag);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Problem Number : %d\n",*EndTag);CHKERRQ(ierr);
      while(*EndTag != 0) {

          ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&alphamat);CHKERRQ(ierr);
          // Set alphamat 
          ierr = MatAXPY(Ap,1,alphamat,SUBSET_NONZERO_PATTERN);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Alphavec is received and set \n");CHKERRQ(ierr);
      
/////////////////////////////// 
/*  Solver Ax =b             */
///////////////////////////////


  ierr = PetscGetTime(&t1);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp1,Ap,Ac,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSolve(ksp1,b,x);CHKERRQ(ierr);
  ierr = PetscGetTime(&t2);CHKERRQ(ierr);
  elapsed_time = t2 - t1;
  ierr = KSPGetIterationNumber(ksp1,&its);CHKERRQ(ierr); 
  PetscPrintf(PETSC_COMM_WORLD,"PETSC: Velocity field is solve in %4.2e secs with %d iterations \n",elapsed_time,its);
  
  ierr = PetscGetTime(&t1);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp2,Ap,Ac,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSolve(ksp2,c,adjvec);CHKERRQ(ierr);
  ierr = PetscGetTime(&t2);CHKERRQ(ierr);
  elapsed_time = t2 - t1;
  ierr = KSPGetIterationNumber(ksp2,&its);CHKERRQ(ierr);  
  PetscPrintf(PETSC_COMM_WORLD,"PETSC: Adjoint equation is solve in %4.2e secs with %d iterations \n",elapsed_time,its);
  
  // Recover Ap to be A
  ierr = MatDuplicate(A,MAT_COPY_VALUES,&Ap);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp1,PETSC_TRUE);CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(ksp2,PETSC_TRUE);CHKERRQ(ierr);
  
 

///////////////////////////////////// 
/*  Find Gradient Vector and fobj  */
/////////////////////////////////////

  // fobj
  ierr = VecDot(c,x,&fobj);CHKERRQ(ierr);
  // gradf
  ierr = VecPointwiseMult(y, x,adjvec);CHKERRQ(ierr);
  ierr = MatMult(Prt,y,gradf);CHKERRQ(ierr);
  
  PetscPrintf(PETSC_COMM_WORLD,"PETSC: fobj anf gradf have been found.\n");

/////////////////////////////// 
/*  Send x to Matlab socket  */
///////////////////////////////

    PetscRealView(1,&fobj, socketviewer);    
    VecView(x,  socketviewer);
    VecView(gradf,  socketviewer); 
 
    PetscPrintf(PETSC_COMM_WORLD,"PETSC: fobj, x, and gradf has been sent to Matlab socket.\n\n");
     
   
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

  ierr = KSPDestroy(ksp1);CHKERRQ(ierr);
  ierr = KSPDestroy(ksp2);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(c);CHKERRQ(ierr);
  ierr = VecDestroy(adjvec);CHKERRQ(ierr);
   
 
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = MatDestroy(Ac);CHKERRQ(ierr);
  ierr = MatDestroy(Ap);CHKERRQ(ierr);  ierr = MatDestroy(alphamat);CHKERRQ(ierr);
  
  

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}



