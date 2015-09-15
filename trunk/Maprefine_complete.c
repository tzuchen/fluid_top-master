static char help[] = "Calculate a Markov martix for Standard map";





#include "petscksp.h"
#include "petscda.h"
#include "petsctime.h"
#include <math.h>
#include <stdlib.h>
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
  PetscInt       rank,size,npt;  
  PetscErrorCode ierr;
  Vec            x,y0,tempvec, *vinda,*vindb,*vindc;
  PetscInt       i,j,k,l,n,p,m,m2,pmax,puse,Istart,Iend,localsize,niter;

  PetscScalar    dx,dy,dx2,dy2;  
  PetscScalar    *Mixnorm;
  PetscInt       iter,*iterind,*nind;
  FILE           *fidoutput;   
  char           fname[50];
  PetscViewer    socketviewer; 
  PetscInt       withMatlab;
  PetscTruth     Matlabflag;


  PetscLogDouble v1,v2,elapsed_time;


  PetscInitialize(&argc,&args,(char *)0,help);
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
  dx2        = dx/2-dx/1e6;
  dy2        = dy/2-dy/1e6;
  npt        = 5;
  pmax       = 4e6;
  puse       = pmax;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: estimated buffer size (per processer) %f Mbytes \n", pmax*1.0/1e6*8*16 );
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: estimated variable size %f Mbytes\n", 1.0*n*n/1e6*8*2);

/////////////////////////////////////////////////////////////////////////////////////
  
  ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE ,n,&tempvec);CHKERRQ(ierr);
  ierr  = VecGetOwnershipRange(tempvec,&Istart,&Iend);CHKERRQ(ierr); 
  localsize = Iend-Istart;
  ierr = VecDestroy(tempvec);CHKERRQ(ierr);
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
    ierr = VecCreateSeq(PETSC_COMM_SELF,pmax*npt,&y0);CHKERRQ(ierr);
     
    ierr =  VecNorm(x0,NORM_2,Mixnorm); CHKERRQ(ierr);  
    PetscPrintf(PETSC_COMM_WORLD,"PETSC: initial norm= %f \n",*(Mixnorm+0)/n ); 

///////////////////////////////////////////////////////////////////////////
// Map Center Points
  PetscInt     *NzindJ,*idx,*idy,*idp;
  PetscScalar  *CenterX,*CenterY,*VecVal,*pty;
  PetscScalar  *ShiftX,*ShiftY,CX,CY, *yarray;
  IS           isx,isy;
  VecScatter   ctx;

  CenterX      = malloc(npt*sizeof(PetscScalar));
  CenterY      = malloc(npt*sizeof(PetscScalar));
  ShiftX       = malloc(npt*sizeof(PetscScalar));
  ShiftY       = malloc(npt*sizeof(PetscScalar));
  VecVal       = malloc(npt*sizeof(PetscScalar));
  yarray       = malloc(pmax*sizeof(PetscScalar)); 

  NzindJ       = malloc(pmax*npt*sizeof(PetscInt));
  idx          = malloc(pmax*npt*sizeof(PetscInt));
  idy          = malloc(pmax*npt*sizeof(PetscInt)); 
  idp          = malloc(pmax*sizeof(PetscInt)); 
 

  *(ShiftX+0) = 0;
  *(ShiftY+0) = 0;
  *(ShiftX+1) = -dx2;
  *(ShiftY+1) = -dy2;
  *(ShiftX+2) =  dx2;
  *(ShiftY+2) = -dy2;
  *(ShiftX+3) = -dx2;
  *(ShiftY+3) =  dy2;
  *(ShiftX+4) =  dy2;
  *(ShiftY+4) =  dx2;


  //*(ShiftX+5) = 0;
  //*(ShiftY+5) = -dy2;
  //*(ShiftX+6) = -dx2;
  //*(ShiftY+6) = 0;
  //*(ShiftX+7) =  dx2;
  //*(ShiftY+7) = 0;
  //*(ShiftX+8) = 0;
  //*(ShiftY+9) =  dy2;

  for(i=0;i<npt*pmax;i++){ *(idy+i)=i; }
  ISCreateGeneralWithArray(PETSC_COMM_SELF,npt*pmax,idy,&isy);
  vinda = &x0;
  vindb = &x;


   sprintf(fname, "mixnorm_%d_%d",n,iter);
   ierr =PetscPrintf(PETSC_COMM_WORLD,"\n iter     norm      time      unit time\n");CHKERRQ(ierr);
   ierr =PetscFOpen(PETSC_COMM_WORLD,fname,"w",&fidoutput);CHKERRQ(ierr);


for(niter=0;niter<iter;niter++){
   
 ierr = PetscGetTime(&v1);CHKERRQ(ierr);
 l = 0; p = 0;

 if (n*localsize-l<=pmax){puse = n*localsize-l;}else{puse=pmax;}     
     for(i=Istart;i<Iend;i++){
         for(j=0;j<n;j++){
              CX = dx2+i*dx;
              CY = dy2+j*dy;              
              for(k=0;k<npt;k++){ 
                   *(CenterX+k) = CX + *(ShiftX+k);
                   *(CenterY+k) = CY + *(ShiftY+k);    
                   InverseStandardMap((CenterX+k),(CenterY+k));
                   *(NzindJ+p*npt +k) =  floor(*(CenterX+k)*n)*n +  floor(*(CenterY+k)*n);      
              }          
              *(idp+p) = Istart*n+ l;
      
             if(p>=puse-1){ 
          
                 ierr =  ISCreateGeneralWithArray(PETSC_COMM_WORLD,npt*puse,NzindJ,&isx);CHKERRQ(ierr);
                 for(m=0;m<npt*puse;m++){ *(idy+m)=m; }
                 ierr =  ISCreateGeneralWithArray(PETSC_COMM_SELF,npt*puse,idy,&isy);CHKERRQ(ierr);
                 ierr =  VecScatterCreate(*vinda,isx,y0,isy,&ctx);CHKERRQ(ierr);
                 ierr =  VecScatterBegin(*vinda,y0,INSERT_VALUES,SCATTER_FORWARD,ctx);CHKERRQ(ierr);
                 ierr =  VecScatterEnd(*vinda,y0,INSERT_VALUES,SCATTER_FORWARD,ctx);CHKERRQ(ierr);
                 ierr =  VecScatterDestroy(ctx);
                 ierr =  VecGetArray(y0,&pty);CHKERRQ(ierr);
              
                 for(m=0;m<puse;m++){
                     for(m2=0;m2<npt;m2++){
                        *(yarray+m) = *(yarray+m)+*(pty+m*npt+m2);
                     }  
                     *(yarray+m) = *(yarray+m)/npt;
                 } 
                 VecRestoreArray(y0,&pty);
                 VecSetValues(*vindb,puse,idp,yarray,INSERT_VALUES);       

 

                 for(m=0;m<pmax;m++){*(yarray+m) = 0; } 
                 p = 0;

                 if (n*localsize-l<=pmax){puse = n*localsize-l-1;}else{puse=pmax;}            
             }else{p++;}

             l++;
         
        }
    }


   VecAssemblyBegin(*vindb);
   VecAssemblyEnd(*vindb);

   vindc = vindb;
   vindb = vinda;
   vinda = vindc;   


   //ierr =  VecCopy(x,x0);CHKERRQ(ierr);
   ierr =  VecNorm(*vinda,NORM_2,Mixnorm+niter); CHKERRQ(ierr); 
   *(Mixnorm+niter) = *(Mixnorm+niter)/n; 

        
   ierr = PetscGetTime(&v2);CHKERRQ(ierr);
   elapsed_time = v2 - v1; 
   PetscPrintf(PETSC_COMM_WORLD,"     %d   %f   %f  %f \n",niter,*(Mixnorm+niter),elapsed_time,elapsed_time/n/n*1e6 );
   PetscFPrintf(PETSC_COMM_WORLD,fidoutput,"    %d   %f   %f  %f\n"
                ,niter,*(Mixnorm+niter),elapsed_time,elapsed_time/n/n*1e6 );
}



 PetscFClose(PETSC_COMM_WORLD,fidoutput);

///////////////////////////////////////////////////////////////////////////

    if(withMatlab==1){
     VecView(x0,socketviewer);
     PetscScalarView(iter,Mixnorm,socketviewer);
    }
 
  free(CenterX);
  free(CenterY);
  free(ShiftX);
  free(ShiftY);
  

  free(x0array);
  free(idx);
  free(idy);
  free(idp);
  free(yarray);

 

  free(NzindJ);

  free(Mixnorm);

 
   ierr = VecDestroy(x0);CHKERRQ(ierr);
   ierr = VecDestroy(x);CHKERRQ(ierr);
   ierr = VecDestroy(y0);CHKERRQ(ierr);
 
  PetscPrintf(PETSC_COMM_WORLD,"Done!");
  
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
  //if(*x==1){*x=0;}
  //if(*y==1){*y=0;}

return 0;
}











