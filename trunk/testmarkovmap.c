
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

  Vec            sol; 
  Vec            dA; 
  Vec6           vecf, cvec;
  
  Mat            A;
  Mat6           Dmat;

  PetscInt       solsize;
  PetscInt       rank,size;  
 
  PetscErrorCode ierr;
  PetscViewer    socketviewer;
  CoordInfo      *cdinfo;
 
  PetscInt       *ny,*nz;  
  PetscScalar *carray[6];
  PetscInt    nnzA;

   


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
 

  // 
  ierr = IntReceive(socketviewer, &ny);CHKERRQ(ierr); 
  ierr = IntReceive(socketviewer, &nz);CHKERRQ(ierr); 
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: ny and nz are received from Matlab \n");CHKERRQ(ierr); 




 solsize   = cdinfo->nofgrid->x + cdinfo->nofgrid->y + cdinfo->nofgrid->z + 24;
 MarkovMap(sol, *ny, *nz, cdinfo, &A, &Dmat, &vecf, carray, &nnzA);

 
// Output of MarkovMap
  SparseView(A,socketviewer);

  PetscScalarView(nnzA,carray[0],socketviewer);
  PetscScalarView(nnzA,carray[1],socketviewer);
  PetscScalarView(nnzA,carray[2],socketviewer);
  PetscScalarView(nnzA,carray[3],socketviewer);
  PetscScalarView(nnzA,carray[4],socketviewer);
  PetscScalarView(nnzA,carray[5],socketviewer);
   
// Input of dMarkovMap
  VecReceive(socketviewer, VECMPI,&(cvec.vec1));
  VecReceive(socketviewer, VECMPI,&(cvec.vec2));
  VecReceive(socketviewer, VECMPI,&(cvec.vec3));
  VecReceive(socketviewer, VECMPI,&(cvec.vec4));
  VecReceive(socketviewer, VECMPI,&(cvec.vec5));
  VecReceive(socketviewer, VECMPI,&(cvec.vec6));


  dMarkovMap(&dA, cvec,  Dmat);
// Output of dMarkovMap
  VecView(dA,socketviewer);

  VecView(vecf.vec1,socketviewer);
  VecView(vecf.vec2,socketviewer);
  VecView(vecf.vec3,socketviewer);
  VecView(vecf.vec4,socketviewer);
  VecView(vecf.vec5,socketviewer);
  VecView(vecf.vec6,socketviewer);




  

//////////////////////////////////////////////////////////////////////////////////////////////  

  ierr = VecDestroy(sol);CHKERRQ(ierr);
  ierr = VecDestroy(dA);CHKERRQ(ierr);

  ierr = VecDestroy(vecf.vec1);CHKERRQ(ierr);
  ierr = VecDestroy(vecf.vec2);CHKERRQ(ierr);
  ierr = VecDestroy(vecf.vec3);CHKERRQ(ierr);
  ierr = VecDestroy(vecf.vec4);CHKERRQ(ierr);
  ierr = VecDestroy(vecf.vec5);CHKERRQ(ierr);
  ierr = VecDestroy(vecf.vec6);CHKERRQ(ierr);

  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = MatDestroy(Dmat.Mat1);CHKERRQ(ierr);
  ierr = MatDestroy(Dmat.Mat2);CHKERRQ(ierr);
  ierr = MatDestroy(Dmat.Mat3);CHKERRQ(ierr);
  ierr = MatDestroy(Dmat.Mat4);CHKERRQ(ierr);
  ierr = MatDestroy(Dmat.Mat5);CHKERRQ(ierr);
  ierr = MatDestroy(Dmat.Mat6);CHKERRQ(ierr);


  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////





