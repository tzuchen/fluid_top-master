/*
 
  This is a MATLAB Mex program which sends a variable from matlab
  to a file that Petsc can read.

  Usage: err = send2file(filename,data);  
 
        Written by Tzu-Chen Liang  10/21/2005

  Since this is called from Matlab it cannot be compiled with C++.

  To compile this file use the command:  make BOBT=g matlabcodes

*/



#include <stdio.h>
#include "petscsys.h"
#include "petsc.h"   
#include "src/sys/viewer/impls/socket/socket.h"
#include "mex.h"
#include "matrix.h"



#define PETSC_MEX_ERROR(a) {fprintf(stdout,"SEND: %s \n",a); return ;}
/*-----------------------------------------------------------------*/
/*                                                                 */
/*-----------------------------------------------------------------*/
#undef __FUNCT__  
#define __FUNCT__ "mexFunction"
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
 int            type,t,i;
 PetscScalar    *r,*Pr;
 PetscInt       m,n,l,*Jc,*Ir,nnz;
 int            fp;
 char           *str;


  if (nrhs != 2) PETSC_MEX_ERROR("Send requires two input arguments.");

  /* get fp */

   str=mxArrayToString(prhs[0]);
   PetscBinaryOpen(str,FILE_MODE_WRITE,&fp);
  
  /* get size */
  m = (PetscInt)mxGetM(prhs[1]);
  n = (PetscInt)mxGetN(prhs[1]);
  // First send the size of the matrix
  PetscBinaryWrite(fp,&m,1,PETSC_INT,PETSC_FALSE);  
  PetscBinaryWrite(fp,&n,1,PETSC_INT,PETSC_FALSE);
  l = m*n;
 
  if(mxIsSparse(prhs[1])){

  //Sparse matrix
       nnz = mxGetNzmax(prhs[1]);
       Jc  = mxGetJc(prhs[1]);
       Ir  = mxGetIr(prhs[1]); 
       Pr  = mxGetPr(prhs[1]);

        PetscBinaryWrite(fp,&nnz,1,PETSC_INT,PETSC_FALSE); 
        PetscBinaryWrite(fp,Jc,n+1,PETSC_INT,PETSC_FALSE); 
        PetscBinaryWrite(fp,Ir,nnz,PETSC_INT,PETSC_FALSE);
        PetscBinaryWrite(fp,Pr,nnz,PETSC_SCALAR,PETSC_FALSE);

}else{ 

   // Dense matrix or a vector
   r = mxGetData(prhs[1]);
   // Then we send the matrix
   PetscBinaryWrite(fp,r,l,PETSC_SCALAR,PETSC_FALSE); 


}
  PetscBinaryClose(fp);

  return;
}





