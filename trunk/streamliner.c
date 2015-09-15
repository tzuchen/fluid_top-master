
/* Program usage:  mpirun -np <procs> fluidsolver [-help] [all PETSc options] */ 

static char help[] = "Calculate streamlines";





#include "petscksp.h"
#include "petscda.h"
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

  Vec            sol,x0,y0,z0; 

  PetscInt       rank,size,Istart,Iend,Lsize;
  Mat            Dy,Dz;
 
  PetscErrorCode ierr;
  PetscViewer    socketviewer;
  CoordInfo      *cdinfo;

  // for streamline calculation
  PetscScalar *xarray, *yarray,*zarray;    
  PetscScalar  *x,*y,*z;
   


   PetscInt n = 5500;

   


  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
  // Rank 0 connects to socket, use default socket
  PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: socket opened! \n");CHKERRQ(ierr);  

  // receive cdinfo from Matlab
  ierr = CdinfoReceive(socketviewer, &cdinfo);

  // receive the velocity field from matlab.
  ierr = VecReceive(socketviewer,VECSEQ,&sol);CHKERRQ(ierr); 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: sol is received from Matlab \n");CHKERRQ(ierr); 
 

  ierr = VecReceive(socketviewer,VECMPI,&x0);CHKERRQ(ierr); 
  ierr = VecReceive(socketviewer,VECMPI,&y0);CHKERRQ(ierr); 
  ierr = VecReceive(socketviewer,VECMPI,&z0);CHKERRQ(ierr); 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: x0,y0 and z0 are received from Matlab \n");CHKERRQ(ierr);   
 

   //StreamlineMap(&x0, &y0, &z0, sol, cdinfo, &Dy, &Dz,32, n ,10,&x,&y,&z);
   StreamlineMapr(&x0, &y0, &z0, sol, cdinfo, &Dy, &Dz,32, n ,10);
   
   ierr = VecGetOwnershipRange(x0,&Istart,&Iend);CHKERRQ(ierr); 
   Lsize = Iend - Istart;

 

   VecView(x0,socketviewer);
   VecView(y0,socketviewer);
   VecView(z0,socketviewer);
 
   //PetscScalarView(Lsize*n,x,socketviewer);
   //PetscScalarView(Lsize*n,y,socketviewer);
   //PetscScalarView(Lsize*n,z,socketviewer);


   SparseView(Dy,socketviewer);
   SparseView(Dz,socketviewer);




   ierr = VecRestoreArray(x0,&xarray);CHKERRQ(ierr);
   ierr = VecRestoreArray(y0,&yarray);CHKERRQ(ierr);
   ierr = VecRestoreArray(z0,&zarray);CHKERRQ(ierr);

   ierr = VecDestroy(x0);CHKERRQ(ierr);
   ierr = VecDestroy(y0);CHKERRQ(ierr);
   ierr = VecDestroy(z0);CHKERRQ(ierr);



   ierr = MatDestroy(Dy);CHKERRQ(ierr);
   ierr = MatDestroy(Dz);CHKERRQ(ierr);

 


  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////





