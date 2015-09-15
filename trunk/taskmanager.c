
/* Program usage:  mpirun -np <procs> fluidsolver [-help] [all PETSc options] */ 

static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
  -random_exact_sol		   	 : use a random exact solution vector\n\
  -view_exact_sol   		   	 : write exact solution vector to stdout\n\
  -rtol <error tolerance>          	 : error tolerance\n\   ";




#include "petscksp.h"
#include "petscda.h"
#include <sys/stat.h>
#include <time.h>
#include "MeshStruc.h"
#include "PetscReceive.h"
#include "PetscStreamline.h"

int cornergen(Vec, Vec, Vec,PetscScalar, PetscScalar,Vec, Vec, Vec);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            b,c,x,y,adjvec,gradf,sol, ix,iy,iz,ixcor,iycor,izcor;
  Vec            jx,jy,jz,jxcor,jycor,jzcor;   			
  Mat            A,Ac,Ap,Prt,alphamat,Dyf,Dzf;        		
  KSP            ksp1,ksp2;    	
  PetscInt       m,n,rank,size, *EndTag,its, *ny,*nz;
  PetscErrorCode ierr;
  PetscScalar    rtol=1.e-5,fobj, solsize,*dy,*dz;  
  PetscViewer    socketviewer;
  PetscLogDouble t1,t2,elapsed_time;
  CoordInfo      *cdinfo;
  
  Mat  	         MarkovA;
  Mat6 	         Dmat;
  Vec6           vecf,cvec;
  PetscInt       nnzA;
  PetscScalar    *carray[6];
  Vec            dA;


  PetscInitialize(&argc,&args,(char *)0,help);
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
  // Rank 0 connects to socket, use default socket
     PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer); 
    //PetscViewerSocketOpen(PETSC_COMM_WORLD,0,5050,&socketviewer); 
     ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: socket opened! \n");CHKERRQ(ierr);     
  
  // List of parameters
     ierr = PetscOptionsGetScalar(PETSC_NULL,"-rtol",&rtol,PETSC_NULL);CHKERRQ(ierr);

  
/////////////////////////////// 
/*  Receive A and b          */
///////////////////////////////

     


  // Receive matrix A
    ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&A);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix A received.\n");CHKERRQ(ierr);  
  // Receive the preconditioner matrix Ac
    ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&Ac);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix Ac received.\n");CHKERRQ(ierr);  
  // Receive the mapping matrix Prt (the transepose of Pr)
    ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&Prt);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix Prt received.\n");CHKERRQ(ierr);  
  // Receive RHS vector b
    ierr = VecReceive(socketviewer,VECMPI,&b);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: vector b received.\n");CHKERRQ(ierr);   
  // Receive RHS vector c 
    ierr = VecReceive(socketviewer,VECMPI,&c);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: vector c received.\n");CHKERRQ(ierr);  

  // Receive ix,iy,iz
    ierr = VecReceive(socketviewer,VECMPI,&ix);CHKERRQ(ierr); 
    ierr = VecReceive(socketviewer,VECMPI,&iy);CHKERRQ(ierr); 
    ierr = VecReceive(socketviewer,VECMPI,&iz);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: ix,iy,iz received.\n");CHKERRQ(ierr); 

  // Receive ixcor,iycor,izcor
    ierr = VecReceive(socketviewer,VECMPI,&ixcor);CHKERRQ(ierr); 
    ierr = VecReceive(socketviewer,VECMPI,&iycor);CHKERRQ(ierr); 
    ierr = VecReceive(socketviewer,VECMPI,&izcor);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: ixcor,iycor,izcor received.\n");CHKERRQ(ierr); 


  // Receive dy,dz
    ierr = ScalarReceive(socketviewer,&dy);CHKERRQ(ierr);  
    ierr = ScalarReceive(socketviewer,&dz);CHKERRQ(ierr);  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: dy,dz received.\n");CHKERRQ(ierr); 

  // receive cdinfo from Matlab
    ierr = CdinfoReceive(socketviewer, &cdinfo);


  // Form x,y and adjvec by copying b       
          ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
          ierr = VecDuplicate(b,&y);CHKERRQ(ierr);
          ierr = VecDuplicate(b,&adjvec);CHKERRQ(ierr);

       
  // Form Ap by copying A    
          ierr = MatDuplicate(A,MAT_COPY_VALUES,&Ap);CHKERRQ(ierr);
  // Set ksp
  	  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp1);CHKERRQ(ierr);
  	  ierr = KSPSetTolerances(ksp1,rtol,1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
  	  ierr = KSPSetFromOptions(ksp1);CHKERRQ(ierr);

  	  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp2);CHKERRQ(ierr);
          ierr = KSPSetTolerances(ksp2,rtol,1.e-50,PETSC_DEFAULT,PETSC_DEFAULT);CHKERRQ(ierr);
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

  // set jx,jy,jz;
  
    ierr  = VecDuplicate(ix,&jx);CHKERRQ(ierr);
    ierr  = VecDuplicate(iy,&jy);CHKERRQ(ierr);
    ierr  = VecDuplicate(iz,&jz);CHKERRQ(ierr);

    ierr  = VecDuplicate(ixcor,&jxcor);CHKERRQ(ierr);
    ierr  = VecDuplicate(iycor,&jycor);CHKERRQ(ierr);
    ierr  = VecDuplicate(izcor,&jzcor);CHKERRQ(ierr);

 /////////////////////////////// 
/*  Receive alphavec         */
/////////////////////////////// 



          ierr = IntReceive(socketviewer,&EndTag);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Task Type : %d\n",*EndTag);CHKERRQ(ierr);
   while(*EndTag != 0) {
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
         if(*EndTag == 1 ){
            ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&alphamat);CHKERRQ(ierr);
            // Set alphamat 
            ierr = MatAXPY(Ap,1,alphamat,SUBSET_NONZERO_PATTERN);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Alphavec is received and set \n");CHKERRQ(ierr);
           
           ierr = PetscGetTime(&t1);CHKERRQ(ierr);
           ierr = KSPSetOperators(ksp1,Ap,Ac,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
           ierr = KSPSolve(ksp1,b,x);CHKERRQ(ierr);
           ierr = PetscGetTime(&t2);CHKERRQ(ierr);
           elapsed_time = t2 - t1;
           ierr = KSPGetIterationNumber(ksp1,&its);CHKERRQ(ierr); 
           PetscPrintf(PETSC_COMM_WORLD,"PETSC: Velocity field is solve in %4.2e secs with %d iterations \n",elapsed_time,its);

             // Recover Ap to be A
           ierr = MatDuplicate(A,MAT_COPY_VALUES,&Ap);CHKERRQ(ierr);
           ierr = KSPSetInitialGuessNonzero(ksp1,PETSC_TRUE);CHKERRQ(ierr);
           ierr = KSPSetInitialGuessNonzero(ksp2,PETSC_TRUE);CHKERRQ(ierr); 
          
           VecView(x,  socketviewer);

        }
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
       if (*EndTag == 2 ){

           ierr = VecReceive(socketviewer,VECSEQ,&sol);CHKERRQ(ierr); 
           ierr = IntReceive(socketviewer, &ny);CHKERRQ(ierr); 
           ierr = IntReceive(socketviewer, &nz);CHKERRQ(ierr); 

           solsize   = cdinfo->nofgrid->x + cdinfo->nofgrid->y + cdinfo->nofgrid->z + 24;
           MarkovMap(sol, *ny, *nz, cdinfo, &MarkovA, &Dmat, &vecf, carray, &nnzA);
 
           // Output of MarkovMap
           SparseView(MarkovA,socketviewer);

           PetscScalarView(nnzA,carray[0],socketviewer);
           PetscScalarView(nnzA,carray[1],socketviewer);
           PetscScalarView(nnzA,carray[2],socketviewer);
           PetscScalarView(nnzA,carray[3],socketviewer);
           PetscScalarView(nnzA,carray[4],socketviewer);
           PetscScalarView(nnzA,carray[5],socketviewer);


          }
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
       if (*EndTag == 3 ){

          VecReceive(socketviewer, VECMPI,&(cvec.vec1));
  	  VecReceive(socketviewer, VECMPI,&(cvec.vec2));
  	  VecReceive(socketviewer, VECMPI,&(cvec.vec3));
  	  VecReceive(socketviewer, VECMPI,&(cvec.vec4));
  	  VecReceive(socketviewer, VECMPI,&(cvec.vec5));
  	  VecReceive(socketviewer, VECMPI,&(cvec.vec6));


          dMarkovMap(&dA, cvec,  Dmat);
          // Output of dMarkovMap
 	  VecView(dA,socketviewer);

          }
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
       if(*EndTag == 4 ){
           ierr = MatReceiveTranspose(socketviewer,MATMPIAIJ,&alphamat);CHKERRQ(ierr);
           ierr = VecReceive(socketviewer,VECMPI,&c);CHKERRQ(ierr);
           ierr = VecReceive(socketviewer,VECMPI,&x);CHKERRQ(ierr);

        
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

           // fobj
           ierr = VecDot(c,x,&fobj);CHKERRQ(ierr);
           // gradf
           ierr = VecPointwiseMult(y, x,adjvec);CHKERRQ(ierr);
           ierr = MatMult(Prt,y,gradf);CHKERRQ(ierr);
  
           PetscPrintf(PETSC_COMM_WORLD,"PETSC: fobj anf gradf have been found.\n");

           PetscRealView(1,&fobj, socketviewer);       
           VecView(gradf,  socketviewer); 

        }
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
       if (*EndTag == 5 ){

           ierr = VecReceive(socketviewer,VECSEQ,&sol);CHKERRQ(ierr); 

           cornergen(ix, iy, iz,0 , 0,jx, jy, jz);
           StreamlineMapr(&jx, &jy,   &jz,   sol, cdinfo, &Dyf, &Dzf,32, n ,10);
           VecView(jy,socketviewer);
           VecView(jz,socketviewer);
           SparseView(Dyf,socketviewer);
           SparseView(Dzf,socketviewer);
           PetscPrintf(PETSC_COMM_WORLD,"PETSC: y,z done.\n");
        
         
           ierr  = VecCopy(ixcor,jxcor);CHKERRQ(ierr);
           ierr  = VecCopy(iycor,jycor);CHKERRQ(ierr);
           ierr  = VecCopy(izcor,jzcor);CHKERRQ(ierr);
           StreamlineMapr(&jxcor, &jycor,   &jzcor,   sol, cdinfo, &Dyf, &Dzf,32, n ,10);
           VecView(jycor,socketviewer);
           VecView(jzcor,socketviewer);
           PetscPrintf(PETSC_COMM_WORLD,"PETSC: ycor,zcor done.\n");

         
          } 

  
   
    // Receive EngTag (if EngTag=0, break)
    ierr = IntReceive(socketviewer,&EndTag);CHKERRQ(ierr);
    if (*EndTag < 1) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Received EngTag. Finalize PETSC. \n");CHKERRQ(ierr);
    }else{
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Problem Number : %d\n",*EndTag);CHKERRQ(ierr);
         }

}
/////////////////////////////// 
/*  Destroy Everything       */
///////////////////////////////

  ierr = KSPDestroy(ksp1);CHKERRQ(ierr);
  ierr = KSPDestroy(ksp2);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(c);CHKERRQ(ierr);
  ierr = VecDestroy(adjvec);CHKERRQ(ierr);
   
 
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = MatDestroy(Ac);CHKERRQ(ierr);
  ierr = MatDestroy(Ap);CHKERRQ(ierr);  ierr = MatDestroy(alphamat);CHKERRQ(ierr);
  
  

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
int cornergen(Vec x, Vec y, Vec z,PetscScalar dy , PetscScalar dz, Vec x1, Vec y1, Vec z1)
{   
    PetscErrorCode ierr;

    ierr  = VecCopy(x,x1);CHKERRQ(ierr);
    ierr  = VecCopy(y,y1);CHKERRQ(ierr);
    ierr  = VecCopy(z,z1);CHKERRQ(ierr);
    ierr  = VecShift(y1, dy);CHKERRQ(ierr);
    ierr  = VecShift(z1, dz);CHKERRQ(ierr);


  return 0;
}



