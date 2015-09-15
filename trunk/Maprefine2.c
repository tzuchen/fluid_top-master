/*
This routine performs large-scale simulations of standardmap



Tzu-Chen Liang 11-6-2006
*/

static char help[] = "Perform iterations";

#include "petscksp.h"
#include "petscda.h"
#include "petsctime.h"
#include <math.h>
#include <stdlib.h>
#include "petscbt.h"
#include "MeshStruc.h"
#include "PetscReceive.h"
#include "PetscStreamline.h"
#include "FFTRoutine.h"
#include "IterationFunction.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  PetscInt       rank,size,npt;  
  PetscScalar    dx,dy,cx,cy;
  PetscErrorCode ierr;
  Vec            x,x0,tempvec, *vinda,*vindb,*vindc;
  PetscInt       i,j,k,n,n2,pmax,puse,Istart,Iend,localsize,niter;
  PetscScalar    **x0array, **aarray,**barray;
  PetscInt       *cacheInt;
  PetscScalar    *cacheScalar;  
  DA             myDA;

  PetscScalar    *Mixnorm;
  PetscInt       iter,*iterind,*nind;
  FILE           *fidoutput, *fidtimelog;   
  char           fname[50],ftimelog[50];
  PetscViewer    socketviewer; 
  PetscInt       withMatlab, doFFT, doSmoothing;
  PetscTruth     Matlabflag, FFTflag, Smoothingflag;
  PetscInt       timelogcount;  
  MPI_Status     status;

  PetscLogDouble v1,v2,elapsed_time;
  timelogcount = 0;
 
 
  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-withMatlab",&withMatlab,&Matlabflag);CHKERRQ(ierr);  
  if (Matlabflag == PETSC_FALSE){withMatlab = 0;}else{withMatlab = 1;}
  ierr = PetscOptionsGetInt(PETSC_NULL,"-doFFT",&doFFT,&FFTflag);CHKERRQ(ierr);  
  if (FFTflag == PETSC_FALSE){doFFT = 0;}else{doFFT = 1;}

  ierr = PetscOptionsGetInt(PETSC_NULL,"-doSmoothing",&doSmoothing,&Smoothingflag);CHKERRQ(ierr);  
  if (Smoothingflag == PETSC_FALSE){doSmoothing = 0;}else{doSmoothing = 1;}

  if(withMatlab==1){
  // Rank 0 connects to socket, use default socket
  PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: socket opened! \n");CHKERRQ(ierr); 

  // Receive n from Matlab
  IntReceive(socketviewer, &nind);
  n = *nind;
  //  Receive iter from Matlab
  IntReceive(socketviewer, &iterind);
  iter = *iterind;
 
  }else{
  ierr = PetscOptionsGetInt(PETSC_NULL,"-ngrid",&n,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-niter",&iter,PETSC_NULL);CHKERRQ(ierr);
  }

  
 
/////////////////////////////////////////////////////////////////////////////////////

  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: number of grid is %d \n", n);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: number of iteration is %d \n", iter);


  Mixnorm    = malloc(iter*sizeof(PetscScalar));

  dx         = 1.0/n;
  dy         = 1.0/n;
  n2         = (PetscInt)(n*0.5);
  npt        = 5;
  pmax       = 5e5;
  puse       = pmax;


  PetscInt      logmax = 1000;
  PetscScalar   Timelog[logmax];
  PetscLogDouble t1,t2;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: estimated buffer size (per processer) %f Mbytes \n", pmax*1.0/1e6*8*17 );
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: estimated variable size %f Mbytes\n", 1.0*n*n/1e6*8*1);

/////////////////////////////////////////////////////////////////////////////////////
  
 // ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE ,n,&tempvec);CHKERRQ(ierr);
 // ierr  = VecGetOwnershipRange(tempvec,&Istart,&Iend);CHKERRQ(ierr); 
 // localsize = Iend-Istart;
 // ierr = VecDestroy(tempvec);CHKERRQ(ierr);


/////////////////////////////////////////////////////////////////////////////////////

if(doSmoothing==1){
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n\n\n\nPETSC: Now Do  DACreate2d \n\n\n\n" );
      ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n\n\n\nPETSC: %d  %d  %d\n\n\n\n",n2,n,size);
      DACreate2d(MPI_COMM_WORLD,DA_XYPERIODIC,DA_STENCIL_BOX,n2,n,1,size,1,2,PETSC_NULL,PETSC_NULL,&myDA);
      DACreateGlobalVector(myDA,&x0);
      DAGetCorners(myDA,PETSC_NULL,&Istart,PETSC_NULL,PETSC_NULL,&localsize,PETSC_NULL);
      Iend = Istart+localsize;
}else{
      ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE ,n,&tempvec);CHKERRQ(ierr);
      ierr  = VecGetOwnershipRange(tempvec,&Istart,&Iend);CHKERRQ(ierr); 
      localsize = Iend-Istart;
      ierr = VecDestroy(tempvec);CHKERRQ(ierr);
      VecCreateMPI(PETSC_COMM_WORLD,localsize*n2,PETSC_DETERMINE ,&x0);
}


//ierr = PetscPrintf(PETSC_COMM_WORLD,"\n\n\n\n\nPETSC: So far so good\n\n\n\n");



VecGetArray2d(x0,n2,localsize,0,0,&x0array);

// Create initial vector

   
   for(j=0;j<localsize;j++){
         for(i=0;i<n2;i++){
             cx = (Istart+j+0.5)*dx;  
             x0array[i][j] = cos(2*M_PI*cx);
        }
    }

    
VecRestoreArray2d(x0,n2,localsize,0,0,&x0array);


    ierr = VecDuplicate(x0,&x);CHKERRQ(ierr); 
    ierr =  VecNorm(x0,NORM_2,Mixnorm); CHKERRQ(ierr);  
    PetscPrintf(PETSC_COMM_WORLD,"PETSC: initial norm= %f \n",*(Mixnorm+0)/n ); 
    vinda = &x0;
    vindb = &x;

    sprintf(fname, "mixnorm_%d_%d",n,iter);
    ierr =PetscPrintf(PETSC_COMM_WORLD,"\n iter     norm      time      unit time\n");CHKERRQ(ierr);
    ierr =PetscFOpen(PETSC_COMM_WORLD,fname,"w",&fidoutput);CHKERRQ(ierr);

///////////////////////////////////////////////////////////////////////////////////////////////////
// Memory allocation for the iteration scheme 

 //   cacheInt    = malloc(1*pmax*sizeof(PetscInt));
 //   cacheScalar = malloc(2*pmax*sizeof(PetscScalar));

   cacheInt    = malloc(2*pmax*sizeof(PetscInt));
   cacheScalar = malloc(2*pmax*sizeof(PetscScalar));
///////////////////////////////////////////////////////////////////////////////////////////////////
// Iteration here!

for(niter=0;niter<iter;niter++){

  ierr = PetscGetTime(&v1);CHKERRQ(ierr);
 // BackwardAverage(vinda, vindb, cacheInt, cacheScalar, n,  npt, pmax, Istart,Iend);
 // BackwardAverageR(vinda, vindb, cacheInt, cacheScalar, n,  npt, pmax, Istart,Iend);
   BackwardAverageRL(vinda, vindb, cacheInt, cacheScalar, n,  npt, pmax, Istart,Iend);
   vindc = vindb;
   vindb = vinda;
   vinda = vindc;   
// if(doSmoothing==1){Smoothing(vinda, vindb,n, myDA, Istart,Iend);}
  ierr = PetscGetTime(&v2);CHKERRQ(ierr);
   //vindc = vindb;
   //vindb = vinda;
   //vinda = vindc;   
  
   ierr =  VecNorm(*vinda,NORM_2,Mixnorm+niter); CHKERRQ(ierr); 
   *(Mixnorm+niter) = *(Mixnorm+niter)/n;       
   elapsed_time = v2 - v1; 
   PetscPrintf(PETSC_COMM_WORLD,"     %d   %f   %f  %f \n",niter,*(Mixnorm+niter),elapsed_time,elapsed_time/n/n*1e6 );
   PetscFPrintf(PETSC_COMM_WORLD,fidoutput,"    %d   %f   %f  %f\n"
                ,niter,*(Mixnorm+niter),elapsed_time,elapsed_time/n/n*1e6 ); 
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//Change oremtation of vector
VecGetArray2d(*vinda,n2,localsize,0,0,&aarray);
VecGetArray2d(*vindb,localsize,n2,0,0,&barray);  for(j=0;j<localsize;j++){
       for(i=0;i<n2;i++){
             barray[j][i] = aarray[i][j];       
        }
    }  
VecRestoreArray2d(*vinda,n2,localsize,0,0,&aarray);
VecRestoreArray2d(*vindb,localsize,n2,0,0,&barray);
   vindc = vindb;
   vindb = vinda;
   vinda = vindc;  

////////////////////////////////////////////////////////////////////////////////////////////////////
// FFT part
if(doFFT==1){FFT2D(*vinda,*vindb, localsize, n, Istart,Iend, iter,doSmoothing);}
////////////////////////////////////////////////////////////////////////////////////////////////////
/*
   if(rank==0){
   sprintf(ftimelog, "timelog_%d_%d",n,iter);
   fidtimelog = fopen(ftimelog,"w");
   for(i=0;i<timelogcount;i++){  fprintf(fidtimelog,"%f ",Timelog[i]); }  
   fprintf(fidtimelog,"\n ");  

   for(j = 1;j<size;j++){
     MPI_Recv(Timelog,timelogcount,MPI_DOUBLE,j,j,PETSC_COMM_WORLD,&status);
      for(i=0;i<timelogcount;i++){  fprintf(fidtimelog,"%f ",Timelog[i]); }
      fprintf(fidtimelog,"\n ");
     }
   fclose(fidtimelog);
   }else{ 
   MPI_Send(Timelog ,timelogcount,MPI_DOUBLE,0,rank,PETSC_COMM_WORLD);
   } 
   PetscFClose(PETSC_COMM_WORLD,fidoutput);
*/ 
///////////////////////////////////////////////////////////////////////////

    if(withMatlab==1){
     VecView(*vinda,socketviewer);
     PetscScalarView(iter,Mixnorm,socketviewer);
    }
 

  // free(x0array); 
   free(Mixnorm);
   free(cacheInt);
   free(cacheScalar);
 
   ierr = VecDestroy(x0);CHKERRQ(ierr);
   ierr = VecDestroy(x);CHKERRQ(ierr);
   
 
  PetscPrintf(PETSC_COMM_WORLD,"Done!");
  
//////////////////////////////////////////////////////////////////////////////////////
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}










