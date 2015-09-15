
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
#include "MarkovMapFunction.h"




#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{

  Vec            sol,x0,y0,z0; 

  PetscInt       rank,size,Istart,Iend,Lsize;
  Mat            Dy,Dz,MarkovA,dAdy,dAdz;
 
  PetscErrorCode ierr;
  PetscViewer    socketviewer;
  CoordInfo      *cdinfo;

  // for streamline calculation
  PetscScalar *xarray, *yarray,*zarray;    
  PetscScalar  *x,*y,*z;
  PetscInt     n, withMatlab; 
  PetscTruth      Matlabflag,flg;
  PetscInt     ny,nz;
  PetscInt     solmaxp;
  PetscScalar  solmax; 

   n          = 5500;
   withMatlab = 0;
   ny         = 2;
   nz         = 2;


  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// Use Matlab as output 
    ierr = PetscOptionsGetInt(PETSC_NULL,"-withMatlab",&withMatlab,&Matlabflag);CHKERRQ(ierr); 
    ierr = PetscOptionsGetInt(PETSC_NULL,"-ngridy",&ny,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(PETSC_NULL,"-ngridz",&nz,&flg);CHKERRQ(ierr);

 
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

          //test
          ierr =  VecLoadfromFile(PETSC_COMM_WORLD,"sol",VECSEQ,&sol);CHKERRQ(ierr);
          ierr = CdinfoLoad(PETSC_COMM_WORLD,"cdinfofiled", &cdinfo);CHKERRQ(ierr);

          ierr = VecReceive(socketviewer,VECMPI,&x0);CHKERRQ(ierr); 
          ierr = VecReceive(socketviewer,VECMPI,&y0);CHKERRQ(ierr); 
          ierr = VecReceive(socketviewer,VECMPI,&z0);CHKERRQ(ierr); 
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: x0,y0 and z0 are received from Matlab \n");CHKERRQ(ierr);   
 
     }else{

           // Load cdinfo
           ierr = CdinfoLoad(PETSC_COMM_WORLD,"cdinfofiled", &cdinfo);CHKERRQ(ierr);
    
           ierr =  VecLoadfromFile(PETSC_COMM_WORLD,"sol",VECSEQ,&sol);CHKERRQ(ierr);
           // Load x0, y0, z0 as initial points for streamlines  
           ierr =  VecLoadfromFile(PETSC_COMM_WORLD,"x0filed",VECMPI,&x0);CHKERRQ(ierr);
           ierr =  VecLoadfromFile(PETSC_COMM_WORLD,"y0filed",VECMPI,&y0);CHKERRQ(ierr);
           ierr =  VecLoadfromFile(PETSC_COMM_WORLD,"z0filed",VECMPI,&z0);CHKERRQ(ierr);

     }

   ierr = VecMax(sol,&solmaxp,&solmax);CHKERRQ(ierr);

   StreamlineMap(x0, y0, z0, sol, cdinfo, &Dy, &Dz,0.01/solmax, n ,10,&x,&y,&z);
   MarkovMap2d(sol, ny, nz, cdinfo,&MarkovA,&dAdy,&dAdz);   

   ierr = VecGetOwnershipRange(x0,&Istart,&Iend);CHKERRQ(ierr); 
   Lsize = Iend - Istart;

 if(withMatlab==1){

  	 VecView(x0,socketviewer);
         VecView(y0,socketviewer);
   	 VecView(z0,socketviewer);
 	 VecSave(x0,"xefilew");
         VecSave(y0,"yefilew");
   	 VecSave(z0,"zefilew");
 
 	 PetscScalarView(Lsize*n,x,socketviewer);
 	 PetscScalarView(Lsize*n,y,socketviewer);
  	 PetscScalarView(Lsize*n,z,socketviewer);


  	 SparseView(Dy,socketviewer);
  	 SparseView(Dz,socketviewer);
  }else{

         SparseSave(MarkovA,"A");

 	 VecSave(x0,"xefilew");
         VecSave(y0,"yefilew");
   	 VecSave(z0,"zefilew");
   
         Vec vx,vy,vz;

        ierr =  VecCreateMPIWithArray(PETSC_COMM_WORLD,Lsize*n,PETSC_DECIDE,x,&vx);CHKERRQ(ierr); 
        ierr =  VecCreateMPIWithArray(PETSC_COMM_WORLD,Lsize*n,PETSC_DECIDE,y,&vy);CHKERRQ(ierr); 
        ierr =  VecCreateMPIWithArray(PETSC_COMM_WORLD,Lsize*n,PETSC_DECIDE,z,&vz);CHKERRQ(ierr); 

 	VecSave(vx,"xeallfilew");
 	VecSave(vy,"yeallfilew");
  	VecSave(vz,"zeallfilew");

	ierr = VecDestroy(vx);CHKERRQ(ierr);
        ierr = VecDestroy(vy);CHKERRQ(ierr);
	ierr = VecDestroy(vz);CHKERRQ(ierr);


   }



   ierr = VecRestoreArray(x0,&xarray);CHKERRQ(ierr);
   ierr = VecRestoreArray(y0,&yarray);CHKERRQ(ierr);
   ierr = VecRestoreArray(z0,&zarray);CHKERRQ(ierr);

   ierr = VecDestroy(x0);CHKERRQ(ierr);
   ierr = VecDestroy(y0);CHKERRQ(ierr);
   ierr = VecDestroy(z0);CHKERRQ(ierr);



   ierr = MatDestroy(Dy);CHKERRQ(ierr);
   ierr = MatDestroy(Dz);CHKERRQ(ierr);
   ierr = MatDestroy(dAdy);CHKERRQ(ierr);
   ierr = MatDestroy(dAdz);CHKERRQ(ierr);
   ierr = MatDestroy(MarkovA);CHKERRQ(ierr);
  


  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
////////////////////////////////////////////////////////////////////////////////////





