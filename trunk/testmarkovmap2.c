
/* Program usage:  mpirun -np <procs> fluidsolver [-help] [all PETSc options] */ 

static char help[] = "Calculate streamlines";





#include "petscksp.h"
#include "petscda.h"
#include "slepceps.h"
#include <math.h>
#include "src/mat/impls/aij/seq/aij.h"
#include "src/inline/spops.h"
#include "src/inline/dot.h"
#include "petscbt.h"
#include "MeshStruc.h"
#include "PetscReceive.h"
#include "PetscStreamline.h"





#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{

  Vec            sol;   
  Mat            A,dAdy,dAdz;


  PetscInt       solsize;
  PetscInt       rank,size;  
 
  PetscErrorCode ierr;
  PetscViewer    socketviewer;
  CoordInfo      *cdinfo;
 
  PetscInt        ny,nz,*nyy,*nzz;  
  PetscInt        nnzA;
  PetscInt        withMatlab;
  PetscTruth      Matlabflag, flg; 
    

  ny = 50;
  nz = 50;
  withMatlab = 0; 

  //PetscInitialize(&argc,&args,(char *)0,help);
  SlepcInitialize(&argc,&args,(char*)0,help);

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-withMatlab",&withMatlab,&Matlabflag);CHKERRQ(ierr);  
  if (Matlabflag == PETSC_FALSE){withMatlab = 0;}else{withMatlab = 1;}
  if(withMatlab==1){

    // Rank 0 connects to socket, use default socket
    PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer);  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: socket opened! \n");CHKERRQ(ierr);  

     // receive cdinfo from Matlab
     ierr = CdinfoReceive(socketviewer, &cdinfo);

     // receive the velocity field from matlab.
     ierr = VecReceive(socketviewer,VECSEQ,&sol);CHKERRQ(ierr); 
     ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: sol is received from Matlab \n");CHKERRQ(ierr); 
 

      // 
      ierr = IntReceive(socketviewer, &nyy);CHKERRQ(ierr); 
      ierr = IntReceive(socketviewer, &nzz);CHKERRQ(ierr); 
      ny = *nyy;
      nz = *nzz;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: ny and nz are received from Matlab \n");CHKERRQ(ierr); 

   }else{
   //without matlab
     ierr = PetscOptionsGetInt(PETSC_NULL,"-ngridy",&ny,&flg);CHKERRQ(ierr);
     ierr = PetscOptionsGetInt(PETSC_NULL,"-ngridz",&nz,&flg);CHKERRQ(ierr);

     ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Load cdinfofiled and solbinfiles  \n");CHKERRQ(ierr);
     // Load cdinfo
     ierr = CdinfoLoad(PETSC_COMM_WORLD,"cdinfofiled", &cdinfo);CHKERRQ(ierr);
     // Load sol
     ierr = VecLoadfromFile(MPI_COMM_WORLD,"solbinfiles", VECSEQ,&sol);CHKERRQ(ierr);

 // VecLoadfromFile(PETSC_COMM_WORLD,"x0filed",VECMPI,&sol);CHKERRQ(ierr);
 // VecSavebin(sol, "x0filed2");
 // VecLoadfromFile(PETSC_COMM_WORLD,"x0filed2",VECMPI,&sol);CHKERRQ(ierr);
   }

 solsize   = cdinfo->nofgrid->x + cdinfo->nofgrid->y + cdinfo->nofgrid->z + 24;
 
 MarkovMap2d(sol, ny, nz, cdinfo,&A,&dAdy,&dAdz);
 //VecView(sol,0);

 if(withMatlab==1){
   // Output of MarkovMap
   SparseView(A,socketviewer);
 }else{
   SparseSave(A,"Anew");
 }
 



  

//////////////////////////////////////////////////////////////////////////////////////////////  

  ierr = VecDestroy(sol);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);

 

  //ierr = PetscFinalize();CHKERRQ(ierr);
    ierr = SlepcFinalize();CHKERRQ(ierr);
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////





