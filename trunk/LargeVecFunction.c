/*
Here are the functions to deal with MPI vec which has the global size
larger then 2^32. 


Tzu-Chen Liang   11-17-2006
*/

#include "petscksp.h"
#include "petscda.h"
#include "petsctime.h"
#include <math.h>
#include <stdlib.h>
#include "petscbt.h"
#include "LargeVecFunction.h"
/////////////////////////////////////////////////////////////////////////

int LargeVecCreate(Vec *x, PetscInt nvec, Vec xvec[]){

 PetscInt       size,rank;
 PetscErrorCode ierr;
 PetscInt       lsize,mylsize;
 PetscInt       i;
 PetscScalar    *ptx;
 PetscInt       pstart[nvec], pend[nvec];


 MPI_Comm_size(PETSC_COMM_WORLD,&size);
 MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


 LargeVecGetOwnershipRange(x,nvec,pstart,pend);

 ierr =  VecGetArray(*x,&ptx);CHKERRQ(ierr);
 ierr =  VecGetLocalSize(*x,&lsize);CHKERRQ(ierr);

 for (i=0;i<nvec;i++){
   if(rank>=pstart[i]&&rank<pend[i]){mylsize = lsize;}else{mylsize=0;}
   ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD,mylsize,PETSC_DETERMINE,ptx,xvec+i);CHKERRQ(ierr);
 }
 

 return 0;
}
/////////////////////////////////////////////////////////////////////////
int ISCreateGeneralWithIJ(MPI_Comm comm,Vec x, Vec xvec[],PetscInt nvec, PetscInt nrow,PetscInt pnum, PetscInt *I, PetscInt *J,IS ISfrom[], IS ISto[]){

 PetscInt       size,rank;
 PetscErrorCode ierr;
 PetscInt       i,k,count;
 PetscInt       pstart[nvec],pend[nvec];
 PetscInt       cstart[nvec],cend[nvec];
 PetscInt       fromarray[pnum],toarray[pnum];
 

 MPI_Comm_size(PETSC_COMM_WORLD,&size);
 MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

 LargeVecGetOwnershipRange(&x,nvec,pstart,pend);
 LargeVecGetColumnOwnershipRange(xvec,nvec, nrow, cstart, cend);

 for(i=0;i<nvec;i++){
     count = 0;
     for(k=0;k<pnum;k++){
          if(*(I+k)>=cstart[i]&&*(I+k)<cend[i]){
          fromarray[count] =  (PetscInt)((*(I+k)-cstart[i])*nrow + *(J+k));
            toarray[count] =  k;
          count++;
          }
     }
    
     ierr =  ISCreateGeneral(PETSC_COMM_SELF,count,fromarray,ISfrom+i);CHKERRQ(ierr);
     ierr =  ISCreateGeneral(PETSC_COMM_SELF,count,toarray,ISto+i);CHKERRQ(ierr);

     count = 0;
 }


 return 0;
}


/////////////////////////////////////////////////////////////////////////
int LargeVecGetOwnershipRange(Vec *x,PetscInt nvec,PetscInt pstart[], PetscInt pend[]){


 PetscInt       size,rank,lsize;
 PetscInt       i,k;
 
 MPI_Comm_size(PETSC_COMM_WORLD,&size);
 MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

 lsize = (PetscInt)floor(size*1.0/nvec); 
 k = 0;
 for(i=0;i<nvec;i++){ 
    pstart[i] = k;
    k = k+lsize;
    pend[i]   = k;
    if(i==nvec-1){pend[i]   = size;}     
 }


 return 0;
}


/////////////////////////////////////////////////////////////////////////
int LargeVecGetColumnOwnershipRange(Vec xvec[],PetscInt nvec,PetscInt nrow,PetscInt cstart[], PetscInt cend[]){
 
  PetscInt i,ncol[nvec];
  PetscInt size;
  
  cstart[0] = 0;

  for(i=0;i<nvec;i++){
     VecGetSize(xvec[i],&size);
     ncol[i] =  (PetscInt)(size*1.0/nrow);
     cend[i]= cstart[i]+ncol[i];
     if(i<nvec-1){cstart[i+1] = cend[i];}     
  }


 return 0;
}
////////////////////////////////////////////////////////////////////////////
int LargeVecScatterCreate(Vec xvec[],IS isx[],Vec y,IS isy[] ,VecScatter ctx[],PetscInt nvec){
   
  PetscErrorCode ierr;
  PetscInt i;

 for(i=0;i<nvec;i++){
     ierr =  VecScatterCreate(xvec[i],isx[i],y,isy[i],ctx+i);CHKERRQ(ierr);
 }


return 0;
}

////////////////////////////////////////////////////////////////////////////
int  LargeVecScatterBeginEnd(Vec xvec[],Vec y,InsertMode addv,ScatterMode mode,VecScatter ctx[],PetscInt nvec){
 
   PetscErrorCode ierr;
   PetscInt i;    
 
   for(i=0;i<nvec;i++){
      ierr =  VecScatterBegin(xvec[i],y,addv,mode,ctx[i]);CHKERRQ(ierr);
      ierr =  VecScatterEnd(xvec[i],y,addv,mode,ctx[i]);CHKERRQ(ierr);
   }
return 0;
}
///////////////////////////////////////////////////////////////////////////
int ISArrayDestroy(IS isx[],PetscInt size){
   PetscErrorCode ierr;
   PetscInt i;    
  for(i=0;i<size;i++){ierr = ISDestroy(isx[i]);CHKERRQ(ierr);}
return 0;
}
///////////////////////////////////////////////////////////////////////////
int VecScatterArrayDestroy(VecScatter ctx[],PetscInt size){
   PetscErrorCode ierr;
   PetscInt i;    
  for(i=0;i<size;i++){ierr = VecScatterDestroy(ctx[i]);CHKERRQ(ierr);}
return 0;
}
///////////////////////////////////////////////////////////////////////////
int VecArrayDestroy(Vec x[],PetscInt size){
   PetscErrorCode ierr;
   PetscInt i;    
  for(i=0;i<size;i++){ierr = VecDestroy(x[i]);CHKERRQ(ierr);}
return 0;


}






 


