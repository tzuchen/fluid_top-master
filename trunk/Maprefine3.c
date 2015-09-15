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
int StandardMap(PetscScalar*,PetscScalar*);
int SkewSymmetricPoint(PetscInt*, PetscInt*, PetscInt);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  PetscInt       rank,size,npt;  
  PetscErrorCode ierr;
  Vec            x,y0,tempvec, *vinda,*vindb,*vindc;


  PetscInt       i,j,k,l,n,n2,p,m,m2,pmax,puse,Istart,Iend,localsize,niter;
  PetscInt       *pi,*pj,*li,*lj,pcount;

  PetscScalar    dx,dy,dx2,dy2,*ptveca,*ptvecb;  
  PetscScalar    *Mixnorm;
  PetscInt       iter,*iterind,*nind;
  FILE           *fidoutput, *fidtimelog;   
  char           fname[50],ftimelog[50];
  PetscViewer    socketviewer; 
  PetscInt       withMatlab;
  PetscTruth     Matlabflag;
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


  if(withMatlab==1){
  // Rank 0 connects to socket, use default socket
//  PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer); 
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
  npt        = 4;
  pmax       = 5e5;
  puse       = pmax;
  n2         = (PetscInt)(n*0.5);

  PetscInt      logmax = 1000;
  PetscScalar   Timelog[logmax];
  PetscLogDouble t1,t2;

PetscPrintf(PETSC_COMM_WORLD,"%d \n" ,6*iter*((PetscInt)(n*1.0/pmax*5)+1)*((PetscInt)(n*1.0/pmax*5)+1)+10 ); 


  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: estimated buffer size (per processer) %f Mbytes \n", pmax*1.0/1e6*8*16 );
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: estimated variable size %f Mbytes\n", 1.0*n*n/1e6*8*1);

/////////////////////////////////////////////////////////////////////////////////////
  
  ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE ,n,&tempvec);CHKERRQ(ierr);
  ierr  = VecGetOwnershipRange(tempvec,&Istart,&Iend);CHKERRQ(ierr); 
  localsize = Iend-Istart;
  ierr = VecDestroy(tempvec);CHKERRQ(ierr);
/////////////////////////////////////////////////////////////////////////////////////
// Create initial vector
    Vec         x0;
    PetscScalar *x0array;
    x0array = malloc((localsize)*n2*sizeof(PetscScalar));
    k = 0;
    for(i=Istart;i<Iend;i++){
        for(j=0;j<n2;j++){
            *(x0array+k) = cos(2*M_PI*(dx/2+i*dx));
            //*(x0array+k) = cos(2*M_PI*(dy/2+j*dy));       
            k++;

        }
    }

  

    ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,n2*localsize,PETSC_DECIDE,x0array,&x0);CHKERRQ(ierr);
    ierr = VecDuplicate(x0,&x);CHKERRQ(ierr); 
    ierr = VecCreateSeq(PETSC_COMM_SELF,pmax*npt,&y0);CHKERRQ(ierr);
     
  //  VecScale (x0,1.0/npt);   

    ierr =  VecNorm(x0,NORM_2,Mixnorm); CHKERRQ(ierr);  
    PetscPrintf(PETSC_COMM_WORLD,"PETSC: initial norm= %f \n",*(Mixnorm+0)/n ); 

   

///////////////////////////////////////////////////////////////////////////
// Map Center Points
  PetscInt     *NzindJ0,*NzindJ1,*NzindJ2,*NzindJ3,*ShiftX,*ShiftY;
  PetscScalar  CenterX,CenterY,*pty;
  PetscScalar  CX,CY;
  PetscScalar  *we0,*we1,*we2,*we3;
  PetscScalar  xd,yd;
  IS           isx,isy;
  VecScatter   ctx;




  ShiftX       = malloc(npt*sizeof(PetscInt));
  ShiftY       = malloc(npt*sizeof(PetscInt));
 
  pi           = malloc(npt*sizeof(PetscInt));
  pj           = malloc(npt*sizeof(PetscInt));
  li           = malloc(npt*sizeof(PetscScalar));
  lj           = malloc(npt*sizeof(PetscScalar));

  NzindJ0      = malloc(pmax*sizeof(PetscInt));
  NzindJ1      = malloc(pmax*sizeof(PetscInt));
  NzindJ2      = malloc(pmax*sizeof(PetscInt));
  NzindJ3      = malloc(pmax*sizeof(PetscInt));
  we0          = malloc(pmax*sizeof(PetscScalar));
  we1          = malloc(pmax*sizeof(PetscScalar));
  we2          = malloc(pmax*sizeof(PetscScalar));
  we3          = malloc(pmax*sizeof(PetscScalar));


  *(ShiftX+0) = 0;
  *(ShiftY+0) = 0;
  *(ShiftX+1) = 1;
  *(ShiftY+1) = 0;
  *(ShiftX+2) = 0;
  *(ShiftY+2) = 1;
  *(ShiftX+3) = 1;
  *(ShiftY+3) = 1;


  vinda = &x0;
  vindb = &x;


   sprintf(fname, "mixnorm_%d_%d",n,iter);
   ierr =PetscPrintf(PETSC_COMM_WORLD,"\n iter     norm      time      unit time\n");CHKERRQ(ierr);
   ierr =PetscFOpen(PETSC_COMM_WORLD,fname,"w",&fidoutput);CHKERRQ(ierr);


 

PetscGetTime(&t1);

for(niter=0;niter<iter;niter++){
   
 ierr = PetscGetTime(&v1);CHKERRQ(ierr);
 l = 0; p = 0; pcount = 0;
 ierr =  VecGetArray(*vinda,&ptveca);CHKERRQ(ierr);

 
 
 if (n2*localsize-l<=pmax){puse = n2*localsize-l;}else{puse=pmax;}  
 
 for(i=Istart;i<Iend;i++){ 
     for(j=0;j<n2;j++){
              CenterX = dx2+i*dx;
              CenterY = dy2+j*dy; 
              StandardMap(&CenterX,&CenterY);
         
              *(pi+0) = floor((CenterX-dx/2)*n);
              *(pj+0) = floor((CenterY-dy/2)*n);
              xd = (CenterX - (*(pi+0)*dx+dx/2))/dx;
              yd = (CenterY - (*(pj+0)*dx+dy/2))/dy;

              if(*(pi+0)==-1) { *(pi+0)=n-1;}
              if(*(pj+0)==-1) { *(pj+0)=n-1;}
            

              if(*(pi+0)<n-1) {*(pi+1) = *(pi+0) + 1; }else{*(pi+1) = n-1;}    
              *(pj+1) = *(pj+0);
              *(pi+2) = *(pi+0);
               if(*(pj+0)<n-1) {*(pj+2) = *(pj+0) + 1; }else{*(pj+2) = n-1;}    
              *(pi+3) = *(pi+1);
              *(pj+3) = *(pj+2) ;


                //fprintf(stdout,"xd = %f   yd = %f\n",xd,yd );

              *(we0+p) = (1-xd)*(1-yd);              
              *(we1+p) = xd*(1-yd);
              *(we2+p) = (1-xd)*yd;
              *(we3+p) = xd*yd;

              for(k=0;k<npt;k++){ if(*(pj+k)>=n2) {SkewSymmetricPoint(pi+k, pj+k, n);} }
                    
              *(NzindJ0+p) =  (PetscInt)(n2*(*(pi+0)) + *(pj+0));  
              *(NzindJ1+p) =  (PetscInt)(n2*(*(pi+1)) + *(pj+1)); 
              *(NzindJ2+p) =  (PetscInt)(n2*(*(pi+2)) + *(pj+2)); 
              *(NzindJ3+p) =  (PetscInt)(n2*(*(pi+3)) + *(pj+3));                                        
              
                          
              if(p>=puse-1){   

                     for(k=0;k<puse;k++){
                        *(we0+k) = *(we0+k)* (*(ptveca+pcount+k));
                        *(we1+k) = *(we1+k)* (*(ptveca+pcount+k));
                        *(we2+k) = *(we2+k)* (*(ptveca+pcount+k));
                        *(we3+k) = *(we3+k)* (*(ptveca+pcount+k));

                      }
                                 
                     VecSetValues(*vindb,puse,NzindJ0,we0,ADD_VALUES);
                     VecSetValues(*vindb,puse,NzindJ1,we1,ADD_VALUES);
                     VecSetValues(*vindb,puse,NzindJ2,we2,ADD_VALUES);
                     VecSetValues(*vindb,puse,NzindJ3,we3,ADD_VALUES);
                  
                     if (n2*localsize-l<=pmax){puse = n2*localsize-l-1;}else{puse=pmax;} 
                     pcount = pcount+p;
                     p = 0;    
                     VecAssemblyBegin(*vindb);
                     VecAssemblyEnd(*vindb);                          
               }else{p++;}
             l++;
       }   
    }
        l=0;
        p=0;
        pcount = 0;

 

 
   VecRestoreArray(*vinda,&ptveca);

   //VecView(*vindb,0);
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
 //  fprintf(fidtimelog,"%d \n", rank);
 
     
   if(niter<iter-1){
    VecScale (*vindb,0); 
    //VecScale (*vinda,1.0/npt);   
    }
}


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
 
///////////////////////////////////////////////////////////////////////////
 
    if(withMatlab==1){
     VecView(x0,socketviewer);
     PetscScalarView(iter,Mixnorm,socketviewer);
    }

  //free(CenterX);
  //free(CenterY);
  free(ShiftX);
  free(ShiftY);
 

  //free(x0array);
 
 

 

  free(NzindJ0);
  free(NzindJ1);
  free(NzindJ2);
  free(NzindJ3);
  free(we0);
  free(we1);
  free(we2);
  free(we3);

  free(Mixnorm);

 
   ierr = VecDestroy(x0);CHKERRQ(ierr);
   ierr = VecDestroy(x);CHKERRQ(ierr);
  // ierr = VecDestroy(y0);CHKERRQ(ierr);
 
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

////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////

int StandardMap(PetscScalar *x,PetscScalar *y){
 
  PetscScalar c;

  c = 0.5;

  *y = fmod(*y+c*sin(2*M_PI*(*x)),1);
  *x = fmod(*x+*y,1);


  if(*x<0){*x=*x+1;}
  if(*y<0){*y=*y+1;}
  //if(*x==1){*x=0;}
  //if(*y==1){*y=0;}

return 0;
}

////////////////////////////////////////////////////////////////////////////

int SkewSymmetricPoint(PetscInt *i, PetscInt *j, PetscInt n){
   *i = n - *i - 1;
   *j = n - *j - 1;
return 0;
}








