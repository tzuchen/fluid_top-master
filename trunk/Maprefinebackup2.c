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
  PetscInt       i,j,k,l,n,*nind,Istart,Iend,localsize;
  PetscScalar    dx,dy,dx2,dy2;  
  PetscScalar    *Mixnorm;
  PetscInt       iter,*iterind;
  PetscInt       npt;  


  
  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
  // Rank 0 connects to socket, use default socket
  PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: socket opened! \n");CHKERRQ(ierr);  
/////////////////////////////////////////////////////////////////////////////////////
 
  // Receive n from Matlab
  IntReceive(socketviewer, &nind);
  n = *nind;
  //  n = 500;
  //  Receive iter from Matlab
  IntReceive(socketviewer, &iterind);
  iter = *iterind;
  //  iter = 2;



  Mixnorm    = malloc(iter*sizeof(PetscScalar));
  dx         = 1.0/n;
  dy         = 1.0/n;
  dx2        = dx/2;
  dy2        = dy/2;
  npt        = 5;

/////////////////////////////////////////////////////////////////////////////////////
  
  ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE ,n,&tempvec);CHKERRQ(ierr);
  ierr  = VecGetOwnershipRange(tempvec,&Istart,&Iend);CHKERRQ(ierr); 
  localsize = Iend-Istart;



  MatCreateMPIAIJ(PETSC_COMM_WORLD,n*localsize,n*localsize,PETSC_DETERMINE,PETSC_DETERMINE,npt,PETSC_NULL,npt,PETSC_NULL,&A);CHKERRQ(ierr);

// Map Center Points
  PetscInt *NzindI, *NzindJ;
  PetscScalar *CenterX,*CenterY, *MapedCenterX,*MapedCenterY,*MatVal;
  

  CenterX      = malloc(9*sizeof(PetscScalar));
  CenterY      = malloc(9*sizeof(PetscScalar));
  MatVal       = malloc(9*sizeof(PetscScalar));
  NzindI       = malloc(5*sizeof(PetscInt));
  NzindJ       = malloc(9*sizeof(PetscInt));
 

  l = 0;
 
   for(i=0;i<npt;i++){*(MatVal+i)=1;}

  //  *(MatVal+0)=0.999;
  //   for(i=1;i<npt;i++){*(MatVal+i)=0.001/4;}
  
  for(i=Istart;i<Iend;i++){
      for(j=0;j<n;j++){
 
            *(CenterX+0) = dx2+i*dx;
            *(CenterY+0) = dy2+j*dy;
            *(CenterX+1) = dx2+i*dx-dx2;
            *(CenterY+1) = dy2+j*dy-dy2;
            *(CenterX+2) = dx2+i*dx+dx2;
            *(CenterY+2) = dy2+j*dy-dy2;
            *(CenterX+3) = dx2+i*dx-dx2;
            *(CenterY+3) = dy2+j*dy+dy2;
            *(CenterX+4) = dx2+i*dx+dx2;
            *(CenterY+4) = dy2+j*dy+dy2; 
           // *(CenterX+5) = dx2+i*dx;
           // *(CenterY+5) = dy2+j*dy-dy2;
           // *(CenterX+6) = dx2+i*dx-dx2;
           // *(CenterY+6) = dy2+j*dy;
           // *(CenterX+7) = dx2+i*dx+dx2;
           // *(CenterY+7) = dy2+j*dy;
           // *(CenterX+8) = dx2+i*dx;
           // *(CenterY+9) = dy2+j*dy+dy2;
   
          for(k=0;k<npt;k++){     
             InverseStandardMap(CenterX+k,CenterY+k);
             *(NzindJ+k) =  floor(*(CenterX+k)*n)*n +  floor(*(CenterY+k)*n);      
             }

            *(NzindI) =  Istart*n+l;
             //MatSetValues(A,1,NzindI,npt,NzindJ,MatVal,ADD_VALUES);
             MatSetValues(A,1,NzindI,npt,NzindJ,MatVal,INSERT_VALUES);
            l++;


        }
  }


 MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
 MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
 





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
    //  *(nnzRow+k) = 1.0;
  }


  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,n*localsize,PETSC_DECIDE,nnzRow,&norvec);CHKERRQ(ierr);
  ierr = MatDiagonalScale(A,norvec,PETSC_NULL);CHKERRQ(ierr);






///////////////////////////////////////////////////////////////////////////


    Vec         x0,y0,*xind,*yind,*zind;
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

  // SparseView(A,0);
  //  VecView(y0,socketviewer);
     PetscScalarView(iter,Mixnorm,socketviewer);
 
  free(CenterX);
  free(CenterY);
  free(x0array);
  free(nnzRow);


  free(NzindI);
  free(NzindJ);
  free(MatVal);
  free(Mixnorm);

   ierr = MatDestroy(A);CHKERRQ(ierr);
   ierr = VecDestroy(x0);CHKERRQ(ierr);
   ierr = VecDestroy(y0);CHKERRQ(ierr);


//////////////////////////////////////////////////////////////////////////////////////
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}


//////////////////////////////////////////////////////////////////////////////////////

int InverseStandardMap(PetscScalar *x,PetscScalar *y){
 
  PetscScalar c;

  c = 0.5;

  *x = fmod(*x-*y,1);
  *y = fmod(*y-c*sin(2*M_PI*(*x)),1);

  if(*x<0){*x=*x+1;}
  if(*y<0){*y=*y+1;}
  if(*x==1){*x=0;}
  if(*y==1){*y=0;}

return 0;
}











