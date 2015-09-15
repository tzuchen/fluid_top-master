
/* Program usage:  mpirun -np <procs> fluidsolver [-help] [all PETSc options] */ 

static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
  -random_exact_sol		   	 : use a random exact solution vector\n\
  -view_exact_sol   		   	 : write exact solution vector to stdout\n\
  -rtol <error tolerance>          	 : error tolerance\n\   ";




#include "petscksp.h"
#include "petscda.h"
#include "slepceps.h"
#include <sys/stat.h>
#include <time.h>
#include <stdlib.h>
#include "MeshStruc.h"
#include "PetscReceive.h"
#include "PetscStreamline.h"
#include "MarkovMapFunction.h"


#undef __FUNCT__  
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            b,c,x,y,adjvec,gradf;   			
  Mat            A,Ac,Ap,Prt,alphamat,R;        		
  KSP            ksp1,ksp2;    	
  PetscInt       m,n,rank,size, *EndTag,its;
  PetscErrorCode ierr;
  PetscScalar    rtol=1.e-5,fobj;  
  PetscViewer    socketviewer;
  PetscLogDouble t1,t2,elapsed_time;
  
  //PetscInitialize(&argc,&args,(char *)0,help);
  SlepcInitialize(&argc,&args,(char*)0,help);

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
      
  // List of parameters
     ierr = PetscOptionsGetScalar(PETSC_NULL,"-rtol",&rtol,PETSC_NULL);CHKERRQ(ierr);

  
/////////////////////////////// 
/*  Load Data                */
///////////////////////////////
/* The following files are required
    Apfiled,Acfiled,Ppfiled,Rfiled
    bpfiled,cpfiled 
    alphafile
    cdinfofiled
    x0filed,y0filed,z0filed
*/

  // Load matrix A
   ierr = MatLoadTransposefromFile(PETSC_COMM_WORLD,"Apfiled",MATMPIAIJ,&A);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix A loaded.\n");CHKERRQ(ierr);  
  // Load the preconditioner matrix Ac
    ierr = MatLoadTransposefromFile(PETSC_COMM_WORLD,"Acfiled",MATMPIAIJ,&Ac);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix Ac received.\n");CHKERRQ(ierr);  
  // Load the mapping matrix Prt (the transepose of Pr)
    ierr = MatLoadTransposefromFile(PETSC_COMM_WORLD,"Ppfiled",MATMPIAIJ,&Prt);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix Prt received.\n");CHKERRQ(ierr);  

  // Load the permutation matrix R 
    ierr = MatLoadTransposefromFile(PETSC_COMM_WORLD,"Rfiled",MATMPIAIJ,&R);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix R received.\n");CHKERRQ(ierr);  


  // Load RHS vector b
    ierr = VecLoadfromFile(PETSC_COMM_WORLD,"bpfiled",VECMPI,&b);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: vector b received.\n");CHKERRQ(ierr);   
  // Load RHS vector c 
    ierr = VecLoadfromFile(PETSC_COMM_WORLD,"cpfiled",VECMPI,&c);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: vector c received.\n");CHKERRQ(ierr);  

  // Form x,y and adjvec by copying b       
          ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
          ierr = VecDuplicate(b,&y);CHKERRQ(ierr);
          ierr = VecDuplicate(b,&adjvec);CHKERRQ(ierr);

       
  // Form Ap by copying A    
          ierr = MatDuplicate(A,MAT_COPY_VALUES,&Ap);CHKERRQ(ierr);
  // Set ksp
  	  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp1);CHKERRQ(ierr);
  	  ierr = KSPSetTolerances(ksp1,rtol,1.e-4,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  	  ierr = KSPSetFromOptions(ksp1);CHKERRQ(ierr);

  	  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp2);CHKERRQ(ierr);
          ierr = KSPSetTolerances(ksp2,rtol,1.e-4,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
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
      
       
      

          ierr = MatLoadTransposefromFile(PETSC_COMM_WORLD,"alphafiled",MATMPIAIJ,&alphamat);CHKERRQ(ierr);
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
  
  PetscPrintf(PETSC_COMM_WORLD,"PETSC: fobj and gradf have been found.\n");

/////////////////////////////// 
/*  Save x to a binary file  */
///////////////////////////////

    // Reverse permutation
    MatMult(R,x,y);
    PetscViewer binv,ascv;
 
    PetscViewerBinaryOpen(PETSC_COMM_WORLD,"solbin",FILE_MODE_WRITE,&binv);
    VecView(y,  binv);
    
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"sol",&ascv);
    VecView(y,  ascv);
 
    PetscPrintf(PETSC_COMM_WORLD,"PETSC: fobj, x, and gradf have been stored.\n\n");


/////////////////////////////// 
/*  Destroy Everything       */
///////////////////////////////
  PetscPrintf(PETSC_COMM_WORLD,"PETSC:Destroy b,x,c.\n\n");
 
  ierr = KSPDestroy(ksp1);CHKERRQ(ierr);
  ierr = KSPDestroy(ksp2);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(c);CHKERRQ(ierr);
  ierr = VecDestroy(adjvec);CHKERRQ(ierr);
   
 
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = MatDestroy(Ac);CHKERRQ(ierr);
  ierr = MatDestroy(Ap);CHKERRQ(ierr);  ierr = MatDestroy(alphamat);CHKERRQ(ierr);
     
///////////////////////////////////////////////////////////////////////////






   
     Vec            v;
     VecScatter     ctx;   
     PetscInt       sizeg, sizel;
     CoordInfo      *cdinfo;
     Vec            x0,y0,z0;
     Mat            Dy,Dz,MarkovA;
     PetscScalar    *xa,*ya,*za;
     PetscInt       nsize = 500,i,k;
     PetscInt       Lsize,Istart,Iend;
     PetscInt       ny,nz;  
     FILE           *outputfd; 
     MPI_Status      status;

    

   // VecCreateSeq(PETSC_COMM_SELF,m,&v);
    VecScatterCreateToAll(y,&ctx,&v);
    VecScatterBegin(y,v,INSERT_VALUES,SCATTER_FORWARD,ctx);  
    VecScatterEnd(y,v,INSERT_VALUES,SCATTER_FORWARD,ctx); 


   CdinfoLoad(PETSC_COMM_WORLD,"cdinfofiled", &cdinfo);
   VecLoadfromFile(PETSC_COMM_WORLD,"x0filed",VECMPI,&x0);
   VecLoadfromFile(PETSC_COMM_WORLD,"y0filed",VECMPI,&y0);
   VecLoadfromFile(PETSC_COMM_WORLD,"z0filed",VECMPI,&z0);


////////////////////////////////////////////////////////////
  //working here
    ny = 50; nz = 50;

    MarkovMap2(v, ny, nz, cdinfo,&MarkovA);
    
    PetscViewer MatlabViewer;
   // PetscViewerMatlabOpen(PETSC_COMM_WORLD,"MarkovA",FILE_MODE_WRITE,&MatlabViewer);
    PetscViewerASCIIOpen(PETSC_COMM_WORLD,"MarkovA",&MatlabViewer);
    SparseSave(MarkovA,MatlabViewer);
    PetscViewerDestroy(MatlabViewer);

//////////////////////////////////////////////////////////////////////////////////
PetscPrintf(PETSC_COMM_WORLD,"PETSC:Calculating Streamline now.\n\n");
  StreamlineMap(x0, y0, z0, v, cdinfo, &Dy, &Dz,4.0, nsize ,10,&xa,&ya,&za);


   ierr = VecGetOwnershipRange(x0,&Istart,&Iend);CHKERRQ(ierr); 
   Lsize = Iend - Istart;




   
   if(rank==0){
          PetscFOpen(PETSC_COMM_WORLD,"output","w",&outputfd);
          for(i=0;i<Lsize*nsize;i++){
                PetscFPrintf(PETSC_COMM_WORLD,outputfd,"%12.8f %12.8f %12.8f\n",*(xa+i),*(ya+i),*(za+i));
          }
          for(k=1;k<size;k++){
  
             MPI_Recv(&Lsize,1,MPI_INT,k,10*k,PETSC_COMM_WORLD,&status);
      	     MPI_Recv(xa,Lsize*nsize,MPI_DOUBLE,k,10*k+1,PETSC_COMM_WORLD,&status);
      	     MPI_Recv(ya,Lsize*nsize,MPI_DOUBLE,k,10*k+2,PETSC_COMM_WORLD,&status);
      	     MPI_Recv(za,Lsize*nsize,MPI_DOUBLE,k,10*k+3,PETSC_COMM_WORLD,&status);

     	     for(i=0;i<Lsize*nsize;i++){

                PetscFPrintf(PETSC_COMM_WORLD,outputfd,"%12.8f %12.8f %12.8f\n",*(xa+i),*(ya+i),*(za+i));
     	     }

           }
          PetscFClose(PETSC_COMM_WORLD,outputfd);
 
   }else{
      
             MPI_Send(&Lsize,1,MPI_INT,0,rank*10,PETSC_COMM_WORLD);
   	     MPI_Send(xa,Lsize*nsize,MPI_DOUBLE,0,rank*10+1,PETSC_COMM_WORLD);
    	     MPI_Send(ya,Lsize*nsize,MPI_DOUBLE,0,rank*10+2,PETSC_COMM_WORLD);
    	     MPI_Send(za,Lsize*nsize,MPI_DOUBLE,0,rank*10+3,PETSC_COMM_WORLD);

   }




PetscPrintf(PETSC_COMM_WORLD,"PETSC:Destroying v,x0,y0,z0.\n\n");
   free(xa);
   free(ya);
   free(za);

   ierr = VecScatterDestroy(ctx);CHKERRQ(ierr);
   ierr = VecDestroy(v);CHKERRQ(ierr);
   ierr = VecDestroy(x0);CHKERRQ(ierr);
   ierr = VecDestroy(z0);CHKERRQ(ierr);
   ierr = VecDestroy(y0);CHKERRQ(ierr);
  



/////////////////////////////// 
/*  Destroy Everything       */
///////////////////////////////

  
  

  //ierr = PetscFinalize();CHKERRQ(ierr); 
  ierr = SlepcFinalize();CHKERRQ(ierr);
  return 0;
}



