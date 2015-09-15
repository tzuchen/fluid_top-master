static char help[] = "Calculate a Markov martix for Standard map";





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


int InverseStandardMap(PetscScalar*,PetscScalar*);


#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  PetscInt       rank,size;  
  PetscErrorCode ierr;
  PetscViewer    socketviewer;
  Vec            x,y0,tempvec;
  PetscInt       i,j,k,l,n,*nind,Istart,Iend,localsize,niter;
  PetscScalar    dx,dy,dx2,dy2;  
  PetscScalar    *Mixnorm;
  PetscInt       iter,*iterind;
    



  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
  // Rank 0 connects to socket, use default socket
  //PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: socket opened! \n");CHKERRQ(ierr);  
/////////////////////////////////////////////////////////////////////////////////////
 
  // Receive n from Matlab
 // IntReceive(socketviewer, &nind);
 // n = *nind;
    n = 100;
  //  Receive iter from Matlab
 // IntReceive(socketviewer, &iterind);
 // iter = *iterind;
    iter = 5;



  Mixnorm    = malloc(iter*sizeof(PetscScalar));
  dx         = 1.0/n;
  dy         = 1.0/n;
  dx2        = dx/2;
  dy2        = dy/2;


/////////////////////////////////////////////////////////////////////////////////////
  
  ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE ,n,&tempvec);CHKERRQ(ierr);
  ierr  = VecGetOwnershipRange(tempvec,&Istart,&Iend);CHKERRQ(ierr); 
  localsize = Iend-Istart;


/////////////////////////////////////////////////////////////////////////////////////
// Create initial vector
    Vec         x0;
    PetscScalar *x0array;

    x0array = malloc((localsize)*n*sizeof(PetscScalar));
    k = 0;
    for(i=Istart;i<Iend;i++){
        for(j=0;j<n;j++){
            *(x0array+k) = cos(2*M_PI*(dx/2+i*dx));
            //*(x0array+k) = cos(2*M_PI*(dy/2+j*dy));       
            k++;

        }
    }

    
    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,n*localsize,PETSC_DECIDE,x0array,&x0);CHKERRQ(ierr);
    ierr = VecDuplicate(x0,&x);CHKERRQ(ierr); 
    ierr = VecCreateSeq(PETSC_COMM_SELF,9,&y0);CHKERRQ(ierr);
     
///////////////////////////////////////////////////////////////////////////
// Map Center Points
  PetscInt     *NzindI, *NzindJ,*idx,*idy;
  PetscScalar  *CenterX,*CenterY, *MapedCenterX,*MapedCenterY,*VecVal;
  PetscScalar  *ShiftX,*ShiftY,CX,CY;
  PetscScalar  rowsum;
  IS           isx,isy;
  VecScatter   ctx;

  CenterX      = malloc(9*sizeof(PetscScalar));
  CenterY      = malloc(9*sizeof(PetscScalar));
  ShiftX       = malloc(9*sizeof(PetscScalar));
  ShiftY       = malloc(9*sizeof(PetscScalar));
  VecVal       = malloc(9*sizeof(PetscScalar));
  NzindI       = malloc(5*sizeof(PetscInt));
  NzindJ       = malloc(9*sizeof(PetscInt));
  idx          = malloc(9*sizeof(PetscInt));
  idy          = malloc(9*sizeof(PetscInt)); 

  for(i=0;i<9;i++){ *(idy+i)=i; }

  *(ShiftX+0) = 0;
  *(ShiftY+0) = 0;
  *(ShiftX+1) = -dx2;
  *(ShiftY+1) = -dy2;
  *(ShiftX+2) = 0;
  *(ShiftY+2) = -dy2;
  *(ShiftX+3) =  dx2;
  *(ShiftY+3) = -dy2;
  *(ShiftX+4) = -dx2;
  *(ShiftY+4) = 0;
  *(ShiftX+5) =  dx2;
  *(ShiftY+5) = 0;
  *(ShiftX+6) = -dx2;
  *(ShiftY+6) =  dy2;
  *(ShiftX+7) = 0;
  *(ShiftY+7) =  dy2;
  *(ShiftX+8) =  dx2;
  *(ShiftY+8) =  dy2;




  ISCreateGeneralWithArray(PETSC_COMM_SELF,9,idy,&isy);

for(niter=0;niter<iter;niter++){
l = 0;
     for(i=Istart;i<Iend;i++){
         for(j=0;j<n;j++){
  
            CX = dx2+i*dx;
            CY = dy2+j*dy;
          for(k=0;k<9;k++){ 
             *(CenterX+k) = CX + *(ShiftX+k);
             *(CenterY+k) = CY + *(ShiftY+k);    
             InverseStandardMap(CenterX+k,CenterY+k);
             *(NzindJ+k) =  floor(*(CenterX+k)*n)*n +  floor(*(CenterY+k)*n);      
             }
          
    

         ierr =  ISCreateGeneralWithArray(PETSC_COMM_WORLD,9,NzindJ,&isx);CHKERRQ(ierr);
         ierr =  VecScatterCreate(x0,isx,y0,isy,&ctx);CHKERRQ(ierr);
         ierr =  VecScatterBegin(x0,y0,INSERT_VALUES,SCATTER_FORWARD,ctx);CHKERRQ(ierr);
         ierr =  VecScatterEnd(x0,y0,INSERT_VALUES,SCATTER_FORWARD,ctx);CHKERRQ(ierr);
         ierr =  VecScatterDestroy(ctx);

         ierr =  VecSum(y0,&rowsum);CHKERRQ(ierr);
          rowsum = rowsum/9; 
         ierr =  VecSetValue(x,Istart*n+l,rowsum,INSERT_VALUES);CHKERRQ(ierr);
          l++;
       
        }
    }

    VecAssemblyBegin(x);
    VecAssemblyEnd(x);

   ierr =  VecCopy(x,x0);CHKERRQ(ierr);
   ierr =  VecNorm(x0,NORM_2,Mixnorm+niter); CHKERRQ(ierr); 
   *(Mixnorm+niter) = *(Mixnorm+niter)/n; 
   PetscPrintf(PETSC_COMM_WORLD,"\n iter= %d\n",niter );
 
}





///////////////////////////////////////////////////////////////////////////

  // SparseView(A,0);
 //   VecView(x0,socketviewer);
 //    PetscScalarView(iter,Mixnorm,socketviewer);

  PetscPrintf(PETSC_COMM_WORLD,"Done!");
 
  free(CenterX);
  free(CenterY);
  free(x0array);


 

  free(NzindI);
  free(NzindJ);

  free(Mixnorm);
 PetscPrintf(PETSC_COMM_WORLD,"Done2!");
 
   ierr = VecDestroy(x0);CHKERRQ(ierr);
 
  PetscPrintf(PETSC_COMM_WORLD,"Done3!");
  //free(x0array);
//////////////////////////////////////////////////////////////////////////////////////
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}


//////////////////////////////////////////////////////////////////////////////////////

int InverseStandardMap(PetscScalar *x,PetscScalar *y){
 
  PetscScalar c;

  c = 0.1;

  *x = fmod(*x-*y,1);
  *y = fmod(*y-c*sin(2*M_PI*(*x)),1);

  if(*x<0){*x=*x+1;}
  if(*y<0){*y=*y+1;}
  if(*x==1){*x=0;}
  if(*y==1){*y=0;}

return 0;
}


