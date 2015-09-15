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
  Mat            A,localA;
  Vec            x,tempvec;
  PetscInt       i,j,k,n,*nind,Istart,Iend,localsize;
  PetscScalar    dx,dy;  
  PetscScalar    *Mixnorm;
  PetscInt       iter,*iterind;
    



  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
  // Rank 0 connects to socket, use default socket
// PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: socket opened! \n");CHKERRQ(ierr);  
/////////////////////////////////////////////////////////////////////////////////////
 
  // Receive n from Matlab
//  IntReceive(socketviewer, &nind);
//  n = *nind;
    n = 4000;
  // Receive iter from Matlab
//  IntReceive(socketviewer, &iterind);
//  iter = *iterind;
    iter = 2;

  Mixnorm    = malloc(iter*sizeof(PetscScalar));
  dx         = 1.0/n;
  dy         = 1.0/n;


/////////////////////////////////////////////////////////////////////////////////////
  
  ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE ,n,&tempvec);CHKERRQ(ierr);
  ierr  = VecGetOwnershipRange(tempvec,&Istart,&Iend);CHKERRQ(ierr); 
  localsize = Iend-Istart;


// Map Center Points
  PetscInt *NzindI, *NzindJ;
  PetscScalar *CenterX,*CenterY, *MapedCenterX,*MapedCenterY,*MatVal;
  

  NzindI       = malloc((localsize)*n*9*sizeof(PetscInt));
  NzindJ       = malloc((localsize)*n*9*sizeof(PetscInt));
  MatVal       = malloc((localsize)*n*9*sizeof(PetscScalar));
  CenterX      = malloc((localsize)*n*sizeof(PetscScalar));
  CenterY      = malloc((localsize)*n*sizeof(PetscScalar));
  MapedCenterX = malloc((localsize)*n*sizeof(PetscScalar));
  MapedCenterY = malloc((localsize)*n*sizeof(PetscScalar));
  

  k = 0;
  for(i=Istart;i<Iend;i++){
        for(j=0;j<n;j++){
            *(CenterX+k) = dx/2+i*dx;
            *(CenterY+k) = dy/2+j*dy;
            k++;
        }
  }

  ierr = PetscMemcpy(MapedCenterX,CenterX,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscMemcpy(MapedCenterY,CenterY,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);



  //Center
  for(k=0;k<localsize*n;k++){ 
     InverseStandardMap(MapedCenterX+k,MapedCenterY+k); 
     *(NzindI+k*9) =  floor(*(MapedCenterX+k)*n);
     *(NzindJ+k*9) =  floor(*(MapedCenterY+k)*n); 
  }

// UpperLeft
  ierr = PetscMemcpy(MapedCenterX,CenterX,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscMemcpy(MapedCenterY,CenterY,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  for(k=0;k<localsize*n;k++){ 
     *(MapedCenterX+k) = *(MapedCenterX+k)-dx/2; 
     *(MapedCenterY+k) = *(MapedCenterY+k)+dy/2;
     InverseStandardMap(MapedCenterX+k,MapedCenterY+k); 
     *(NzindI+k*9+1) =  floor(*(MapedCenterX+k)*n);
     *(NzindJ+k*9+1) =  floor(*(MapedCenterY+k)*n); 
  }
// UpperCenter
  ierr = PetscMemcpy(MapedCenterX,CenterX,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscMemcpy(MapedCenterY,CenterY,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  for(k=0;k<localsize*n;k++){ 
     
     *(MapedCenterY+k) = *(MapedCenterY+k)+dy/2;
     InverseStandardMap(MapedCenterX+k,MapedCenterY+k); 
     *(NzindI+k*9+2) =  floor(*(MapedCenterX+k)*n);
     *(NzindJ+k*9+2) =  floor(*(MapedCenterY+k)*n); 
  }

// UpperRight
  ierr = PetscMemcpy(MapedCenterX,CenterX,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscMemcpy(MapedCenterY,CenterY,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  for(k=0;k<localsize*n;k++){ 
     *(MapedCenterX+k) = *(MapedCenterX+k)+dx/2; 
     *(MapedCenterY+k) = *(MapedCenterY+k)+dy/2;
     InverseStandardMap(MapedCenterX+k,MapedCenterY+k); 
     *(NzindI+k*9+3) =  floor(*(MapedCenterX+k)*n);
     *(NzindJ+k*9+3) =  floor(*(MapedCenterY+k)*n); 
  }

// CenterLeft
  ierr = PetscMemcpy(MapedCenterX,CenterX,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscMemcpy(MapedCenterY,CenterY,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  for(k=0;k<localsize*n;k++){ 
     *(MapedCenterX+k) = *(MapedCenterX+k)-dx/2; 
     
     InverseStandardMap(MapedCenterX+k,MapedCenterY+k); 
     *(NzindI+k*9+4) =  floor(*(MapedCenterX+k)*n);
     *(NzindJ+k*9+4) =  floor(*(MapedCenterY+k)*n); 
  }

// CenterRight
  ierr = PetscMemcpy(MapedCenterX,CenterX,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscMemcpy(MapedCenterY,CenterY,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  for(k=0;k<localsize*n;k++){ 
     *(MapedCenterX+k) = *(MapedCenterX+k)+dx/2; 
     
     InverseStandardMap(MapedCenterX+k,MapedCenterY+k); 
     *(NzindI+k*9+5) =  floor(*(MapedCenterX+k)*n);
     *(NzindJ+k*9+5) =  floor(*(MapedCenterY+k)*n); 
  }

// DownLeft
  ierr = PetscMemcpy(MapedCenterX,CenterX,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscMemcpy(MapedCenterY,CenterY,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  for(k=0;k<localsize*n;k++){ 
     *(MapedCenterX+k) = *(MapedCenterX+k)-dx/2; 
     *(MapedCenterY+k) = *(MapedCenterY+k)-dy/2;
     InverseStandardMap(MapedCenterX+k,MapedCenterY+k); 
     *(NzindI+k*9+6) =  floor(*(MapedCenterX+k)*n);
     *(NzindJ+k*9+6) =  floor(*(MapedCenterY+k)*n); 
  }

// DownCenter
  ierr = PetscMemcpy(MapedCenterX,CenterX,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscMemcpy(MapedCenterY,CenterY,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  for(k=0;k<localsize*n;k++){ 
     
     *(MapedCenterY+k) = *(MapedCenterY+k)-dy/2;
     InverseStandardMap(MapedCenterX+k,MapedCenterY+k); 
     *(NzindI+k*9+7) =  floor(*(MapedCenterX+k)*n);
     *(NzindJ+k*9+7) =  floor(*(MapedCenterY+k)*n); 
  }

// DownRight
  ierr = PetscMemcpy(MapedCenterX,CenterX,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  ierr = PetscMemcpy(MapedCenterY,CenterY,localsize*n*sizeof(PetscScalar));CHKERRQ(ierr);
  for(k=0;k<localsize*n;k++){ 
     *(MapedCenterX+k) = *(MapedCenterX+k)+dx/2; 
     *(MapedCenterY+k) = *(MapedCenterY+k)-dy/2;
     InverseStandardMap(MapedCenterX+k,MapedCenterY+k); 
     *(NzindI+k*9+8) =  floor(*(MapedCenterX+k)*n);
     *(NzindJ+k*9+8) =  floor(*(MapedCenterY+k)*n); 
  }

///////////////////////////////////////////////////////////////////////////

 for(k=0;k<localsize*n*9;k++){
      *(NzindJ+k) = *(NzindI+k)*n + *(NzindJ+k); 
     // *(NzindI+k) = k/9+Istart*n; 
      *(MatVal+k)= 1;    
 }


///////////////////////////////////////////////////////////////////////////
//Build Matrix A
    PetscInt *Rowind;
    Rowind  = malloc((localsize)*n*sizeof(PetscInt));
    for(k=0;k<(localsize)*n+1;k++){
        *(Rowind+k) = 9*k+0*Istart*n;
    }

 
PetscPrintf(PETSC_COMM_WORLD,"\nBefore create SEQ matrices\n");
      ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,n*localsize,n*n,Rowind,NzindJ,MatVal,&localA);CHKERRQ(ierr); 
PetscPrintf(PETSC_COMM_WORLD,"\nDone creating SEQ matrices\n");


//MatSetValues(Mat mat,PetscInt m,const PetscInt idxm[],PetscInt n,const PetscInt idxn[],const PetscScalar v[],InsertMode addv)
     //MatCreate(PETSC_COMM_WORLD,A);
     //MatSetSizes(A,localsize*n,n*n,PETSC_DETERMINE,PETSC_DETERMINE);
     //MatSetType(A,MATMPIAIJ);
 
    //MatCreateMPIAIJ(PETSC_COMM_WORLD,n*localsize,n*n,PETSC_DETERMINE,PETSC_DETERMINE,9,PETSC_NULL,9,PETSC_NULL,&A);

PetscPrintf(PETSC_COMM_WORLD,"\nDone allocating space for MPI matrix\n");

     //MatMerge(PETSC_COMM_WORLD,localA,PETSC_DECIDE,MAT_REUSE_MATRIX,&A);
      MatMerge(PETSC_COMM_WORLD,localA,PETSC_DECIDE,MAT_INITIAL_MATRIX,&A);

PetscPrintf(PETSC_COMM_WORLD,"\nDoneconverting\n");
///////////////////////////////////////////////////////////////////////////
// Normalize rows
  PetscInt    *Ir,*Jc;
  PetscScalar *Pr,*nnzRow;
  PetscInt    nnz;
  Vec         norvec;
 

  MatGetDataLocal(A,&Ir,&Jc,&Pr,&nnz );
  nnzRow    = malloc((localsize)*n*sizeof(PetscScalar));


  for(k=0;k<(localsize)*n;k++){
    *(nnzRow+k) = (PetscScalar)(1.0/(*(Ir+k+1)-*(Ir+k)));
  }


  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,n*localsize,PETSC_DECIDE,nnzRow,&norvec);CHKERRQ(ierr);
  ierr = MatDiagonalScale(A,norvec,PETSC_NULL);CHKERRQ(ierr);


///////////////////////////////////////////////////////////////////////////

    Vec x0,y0,*xind,*yind,*zind;
    for(k=0;k<(localsize)*n;k++){
    *(CenterX+k) = cos(2*M_PI*(*(CenterX+k)));
    }

    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,n*localsize,PETSC_DECIDE,CenterX,&x0);CHKERRQ(ierr);
    ierr = VecDuplicate(x0,&y0);CHKERRQ(ierr); 
///////////////////////////////////////////////////////////////////////////
    xind = &x0;
    yind = &y0;

   for(k=0;k<iter;k++){
      MatMult(A,*xind,*yind);
      VecNorm(*yind,NORM_2,Mixnorm+k);  
      *(Mixnorm+k) = *(Mixnorm+k)/n; 
      // switch reference
      zind = xind;
      xind = yind;
      yind = zind;  
   }


///////////////////////////////////////////////////////////////////////////

   //SparseView(A,socketviewer);
   //  VecView(y0,socketviewer);
  //  PetscScalarView(iter,Mixnorm,socketviewer);

  PetscPrintf(PETSC_COMM_WORLD,"Done!");
 
  free(CenterX);
  free(CenterY);
  free(MapedCenterX);
  free(MapedCenterY);

  free(nnzRow);
  free(Rowind);

  free(NzindI);
  free(NzindJ);
  free(MatVal);

 
   ierr = MatDestroy(A);CHKERRQ(ierr);

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















