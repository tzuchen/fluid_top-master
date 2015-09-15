/*
This routine performs large-scale simulations of standardmap

  -ngrid           number of grids in each direction
  -niter           number of iterations
  -doFFT           = 1, do FFT after the last iteration
  -doSmoothing     = k, do k times smoothing step after each iteration     
  -doFFTdiffusion  = d, do FFT diffusion, diffusion rate = d
  -Fnameadd        add some string to the rcorded filename

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
#include "LargeVecFunction.h"


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  PetscInt       rank,size,npt;  
  PetscScalar    dx,dy,cx,cy;
  PetscErrorCode ierr;
  Vec            x,x0,tempvec, *vinda,*vindb,*vindc;
  PetscInt       i,j,k,n,n2,pmax,puse,Istart,Iend,localsize,niter;
  PetscScalar    **x0array;
  PetscInt       *cacheInt,*scacheInt;
  PetscScalar    *cacheScalar,*scacheScalar;  
  PetscScalar    vsum;  

  PetscScalar    *Mixnorm1,*Mixnorm2;
  PetscInt       iter,*iterind,*nind;
  FILE           *fidoutput;   
  char           fname[50],fnameadd[50];
  PetscViewer    socketviewer; 
  PetscInt       withMatlab, doFFT, doSmoothing, SmoothingCount=0,IC;
  PetscScalar    doFFTdiffusion,diff;
  PetscTruth     Matlabflag, FFTflag, Smoothingflag, FFTdiffusionflag,Fnameaddflag,ICflag,cflag;
  VecScatter     ctx;
  PetscScalar    c; //c is the exptra parameter to pass to the map, default =0.3

  PetscLogDouble v1,v2,elapsed_time;
  
 
 
  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-withMatlab",&withMatlab,&Matlabflag);CHKERRQ(ierr);  
  if (Matlabflag == PETSC_FALSE){withMatlab = 0;}else{withMatlab = 1;}

  ierr = PetscOptionsGetInt(PETSC_NULL,"-doFFT",&doFFT,&FFTflag);CHKERRQ(ierr);  
  if (FFTflag == PETSC_FALSE){doFFT = 0;}else{/*doFFT = doFFT;*/}

  ierr = PetscOptionsGetInt(PETSC_NULL,"-doSmoothing",&doSmoothing,&Smoothingflag);CHKERRQ(ierr);  
  if (Smoothingflag == PETSC_FALSE){SmoothingCount=0;}else{SmoothingCount = doSmoothing;}

  //diff  = doFFTdiffusion 
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-doFFTdiffusion",&doFFTdiffusion,&FFTdiffusionflag);CHKERRQ(ierr);  
  if (FFTdiffusionflag == PETSC_FALSE){diff = 0;}else{diff = doFFTdiffusion;}

  // Check added file name
  ierr = PetscOptionsGetString(PETSC_NULL,"-Fnameadd",fnameadd,50,&Fnameaddflag);CHKERRQ(ierr); 
  if (Fnameaddflag == PETSC_FALSE){sprintf(fnameadd,"");}

  // Initial Condition: 0=>cos function, 1=>delta function
  ierr = PetscOptionsGetInt(PETSC_NULL,"-IC",&IC,&ICflag);CHKERRQ(ierr);  
  if (ICflag == PETSC_FALSE){IC=0;}

   //diff  = doFFTdiffusion 
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-epsilon",&c,&cflag);CHKERRQ(ierr);  
  if (cflag == PETSC_FALSE){c = 0.3;}else{/* c=c*/}


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


  Mixnorm1    = malloc(iter*sizeof(PetscScalar));
  Mixnorm2    = malloc(iter*sizeof(PetscScalar));

  dx         = 1.0/n;
  dy         = 1.0/n;
  n2         = (PetscInt)(n*0.5);
  npt        = 4;
  pmax       = 5e5;
  puse       = pmax;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: estimated buffer size (per processer) %f Mbytes \n", pmax*1.0/1e6*8*4 );
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: estimated variable size %f Mbytes\n", 1.0*n*n/1e6*8*1);


      ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE ,n,&tempvec);CHKERRQ(ierr);
      ierr  = VecGetOwnershipRange(tempvec,&Istart,&Iend);CHKERRQ(ierr); 
      localsize = Iend-Istart;
      ierr = VecDestroy(tempvec);CHKERRQ(ierr);
      VecCreateMPI(PETSC_COMM_WORLD,localsize*n2,PETSC_DETERMINE ,&x0);


// Create initial vector
   VecGetArray2d(x0,localsize,n2,0,0,&x0array);   
   if((IC==0)||(IC==1)){
   //cos ic
       for(i=0;i<localsize;i++){
            for(j=0;j<n2;j++){
               cx = (Istart+i+0.5)*dx;
               cy = (j+0.5)*dy;
               //if(cx<0.5){x0array[i][j] =-1;}else{x0array[i][j] =1;}       
               if(IC==0){x0array[i][j] = 1+cos(2*M_PI*cx);}
                    else{x0array[i][j] = 1+cos(2*M_PI*cy);}                
            }
       }    }else{ //IC=1
       for(i=0;i<localsize;i++){
            for(j=0;j<n2;j++){     
               x0array[i][j] =0;                
            }
       }
       if(Istart==0){
               for(i=0;i<5;i++){
                  for(j=0;j<5;j++){     
                    x0array[i][j] =1;                
                  }
               }
       }
    }


    VecRestoreArray2d(x0,localsize,n2,0,0,&x0array);

    //normalize 
    ierr =  VecSum(x0,&vsum); CHKERRQ(ierr);
    ierr =  VecScale(x0, n*n*0.5/vsum); CHKERRQ(ierr);

 
    ierr = VecDuplicate(x0,&x);CHKERRQ(ierr); 
    ierr =  VecNorm(x0,NORM_2,Mixnorm2); CHKERRQ(ierr);  
    PetscPrintf(PETSC_COMM_WORLD,"PETSC: initial 2-norm= %f \n",*(Mixnorm2+0) ); 
    vinda = &x0;
    vindb = &x;

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

    if(Smoothingflag==PETSC_FALSE) {sprintf(fname, "mixnorm_%d_%d%s",n,iter,fnameadd);}
                 else  {sprintf(fname, "mixnorm_%d_%d%ss",n,iter,fnameadd);}
 

    ierr =PetscPrintf(PETSC_COMM_WORLD,"\n iter     norm      time      unit time\n");CHKERRQ(ierr);
    ierr =PetscFOpen(PETSC_COMM_WORLD,fname,"w",&fidoutput);CHKERRQ(ierr);

///////////////////////////////////////////////////////////////////////////////////////////////////
// Memory allocation for the iteration scheme 


 //  For BackwardAverage
    cacheInt    = malloc(2*pmax*sizeof(PetscInt));
    cacheScalar = malloc(2*pmax*sizeof(PetscScalar));

 if(Smoothingflag==PETSC_TRUE){ 
      PetscInt bsize; 
      if(n2>localsize){bsize=n2;}else{bsize=localsize;}
      scacheScalar = malloc(4*bsize*sizeof(PetscScalar));
      scacheInt = malloc(8*bsize*sizeof(PetscInt));
     // SkewSymmetricScatter(vinda,scacheScalar ,scacheInt, n2,Istart,localsize ,&ctx);
   }
  

///////////////////////////////////////////////////////////////////////////////////////////////////
// Iteration here!

for(niter=0;niter<iter;niter++){

  //test
  if(FFTdiffusionflag==PETSC_TRUE){FFT2D(*vinda,*vindb, n2, n, Istart,Iend, niter, Smoothingflag,diff,fnameadd);}


   ierr = PetscGetTime(&v1);CHKERRQ(ierr);
   BackwardAverage(vinda, vindb, cacheInt, cacheScalar, n,  npt, pmax, Istart,Iend,c);
   //BackwardCenter(vinda, vindb, cacheInt, cacheScalar, n,  npt, pmax, Istart,Iend);
   vindc = vindb;
   vindb = vinda;
   vinda = vindc;  

 if(Smoothingflag==1){  
     
     for(i=0;i<SmoothingCount;i++){
      SmoothingRL(vinda, vindb,scacheScalar,scacheInt,&ctx,n, Istart,Iend);}
     }
 
    ierr = PetscGetTime(&v2);CHKERRQ(ierr);

// Calculate 1-norm and 2-norm

    ierr =  VecShift(*vinda, -1) ;CHKERRQ(ierr);
    ierr =  VecNorm(*vinda,NORM_1,Mixnorm1+niter); CHKERRQ(ierr);
    *(Mixnorm1+niter) = *(Mixnorm1+niter)/n/n;
    ierr =   VecNorm(*vinda,NORM_2,Mixnorm2+niter); CHKERRQ(ierr); 
    *(Mixnorm2+niter) = *(Mixnorm2+niter)/n;
    ierr =  VecShift(*vinda, 1) ;CHKERRQ(ierr);



   
   elapsed_time = v2 - v1; 
   PetscPrintf(PETSC_COMM_WORLD,"     %d   %20.18f %20.18f  %f  %f \n",niter,*(Mixnorm1+niter),*(Mixnorm2+niter),elapsed_time,elapsed_time/n/n*1e6 );
   PetscFPrintf(PETSC_COMM_WORLD,fidoutput,"    %d   %20.18f %20.18f  %f  %f\n"
                ,niter,*(Mixnorm1+niter),*(Mixnorm2+niter),elapsed_time,elapsed_time/n/n*1e6 ); 



}

//

////////////////////////////////////////////////////////////////////////////////////////////////////
// FFT part
if(doFFT==1){FFT2D(*vinda,*vindb, n2, n, Istart,Iend, iter, Smoothingflag,0,fnameadd);}
////////////////////////////////////////////////////////////////////////////////////////////////////


    if(withMatlab==1){
     VecView(*vinda,socketviewer);
     PetscScalarView(iter,Mixnorm1,socketviewer);
     PetscScalarView(iter,Mixnorm2,socketviewer);
    }
 

 
   free(Mixnorm1);
   free(Mixnorm2);
   free(cacheInt);
   free(cacheScalar);

   if(doSmoothing==1){
   free(scacheScalar);
   free(scacheInt);

   }

   // free(x0array); 
   ierr = VecDestroy(x0);CHKERRQ(ierr);
   ierr = VecDestroy(x);CHKERRQ(ierr);
  
 
  PetscPrintf(PETSC_COMM_WORLD,"Done!");
  
//////////////////////////////////////////////////////////////////////////////////////
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}










