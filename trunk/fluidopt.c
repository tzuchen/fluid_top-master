
/* Program usage:  mpirun -np <procs> fluidopt [-help] [all PETSc options] */ 

static char help[] = "Solves a linear system in parallel with KSP.\n\
Input parameters include:\n\
  -app 1 or 2 or 3            	     : 1 given c \n\  
                                       2 optimize on map \n\
                                       3 minimize the second largest eigenvalue\n\
  -rtol                              : rtol in ksp  \n\
  -abstol                            : abstol in ksp   \n\
  -ngrid                             : number of grids on Markov matrix\n\
  -niter                             : number of iteraiton in optimization \n\
  -stepsize                          : step size \n\
  -operator                          : 0 P-F, 1 Koopman";


#include "petscksp.h"
#include "petscda.h"
#include "slepceps.h"
#include <sys/stat.h>
#include <time.h>
#include <stdlib.h>
#include "MeshStruc.h"
#include "PetscReceive.h"
#include "PetscStreamline.h"
#include "MarkovMapFunction.h"


#undef __FUNCT__  
#define __FUNCT__ "main"
int main(int argc,char **args)
{



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Variables: general
     PetscErrorCode  ierr;
     PetscLogDouble  t1,t2,elapsed_time;
     PetscTruth      flg;
     PetscViewer     binv,ascv;
     PetscInt        withMatlab;
     PetscTruth      Matlabflag;
     PetscViewer     socketviewer;
     PetscInt        iter;
     PetscInt        app,region;
     
// Variables:V solve
     Vec            b,c,cp,x,y,adjvec,gradf,beinp,beinA,regionvec;   			
     Mat            A,Ac,Ap,Prt,alphamat,R;        		
     KSP            ksp1,ksp2;    	
     PetscInt       m,n,rank,size,its,kspmaxits;    
     PetscScalar    rtol,abstol,fobj;  
     PetscScalar    gradfnorm;
     PetscScalar    *gradfarray;
    

// Variables: Markov Matrix
     Vec            v;
     VecScatter     ctx;   
     CoordInfo      *cdinfo;
     Vec            x0,y0,z0;
     Mat            MarkovA,dAdy,dAdz,dydv,dzdv;
     PetscScalar    *xa,*ya,*za;
     PetscInt       nsize,i,j,k;
     PetscInt       Lsize,Istart,Iend;
     PetscInt       ny,nz;  
     FILE           *outputfd; 
     MPI_Status      status;
     PetscInt       operator;

// Variables: Eigenvalue Problem
     EPS            eps;  /* eigenproblem solver context */
     EPSType        type;
     EPSProblemType problemtype;
     PetscReal      epstol;
     PetscScalar    kr, ki;
     Vec            xr, xi;
     int            nev, ncv, maxit, nconv;
     PetscInt       eig;

// Variables: Gradients
     Vec            Ones;
     Vec            dldy, dldz;
     Vec            dldv1,dldv2;
     IS             is1,is2;
     PetscInt       dldvlength;
     VecScatter     ctx1;
     PetscScalar    *beinparray, xrnorm, gradnorm;

//   Output files
     PetscViewer solfileviewer,markovfileviewer;
     char        alphafilename[50];

// Variables: Optimization
     PetscInt       niter;             
     PetscScalar    stepsize; 
     PetscInt       restart;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//Parameters: General
   app     = 1; 
   region  = 0;

//Parameters: Vsolve
   rtol    = 1e-10;
   abstol  = 1e-10;
   kspmaxits = 1e5;
//Parameters: Markov Matrix
   nsize   = 2000;
   ny      = 100;
   nz      = 100;
   eig     = 1;
   operator= 0; //default: P-F operator
   
//parameters: Eigenvalue Problem
   epstol  = 1e-10;
   nev     = 5; 
   ncv     = 30;
   maxit   = 1e6;
   type    = EPSKRYLOVSCHUR;
   problemtype = EPS_NHEP;

//parameters: Optimization
  niter    = 5; 
  stepsize = 5e5;
  restart  = 0;
//////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////// 
//Initialize PETSC 
  //PetscInitialize(&argc,&args,(char *)0,help);
  SlepcInitialize(&argc,&args,(char*)0,help);

  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Petsc Initializes successfully! \n");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: comm_size is %d \n", size);
 
      
  // Get parameters from input
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-rtol",&rtol,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-abstol",&abstol,&flg);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-ngridy",&ny,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-ngridz",&nz,&flg);CHKERRQ(ierr);
 
  ierr = PetscOptionsGetInt(PETSC_NULL,"-niter",&niter,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-app",&app,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetScalar(PETSC_NULL,"-stepsize",&stepsize,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-operator",&operator,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-restart",&restart,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-kspmaxits",&kspmaxits,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-eig",&eig,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-iter",&iter,&flg);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-region",&region,&flg);CHKERRQ(ierr);
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// Use Matlab as output 
    ierr = PetscOptionsGetInt(PETSC_NULL,"-withMatlab",&withMatlab,&Matlabflag);CHKERRQ(ierr);  
    if (Matlabflag == PETSC_FALSE){withMatlab = 0;}else{withMatlab = 1;}
    if(withMatlab==1){
      // Rank 0 connects to socket, use default socket
      PetscViewerSocketOpen(PETSC_COMM_WORLD,0,PETSC_DEFAULT,&socketviewer);  
     ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Matlab socket opened! \n");CHKERRQ(ierr); 
    }

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//Load Data                
/* The following files are required
    Apfiled,Acfiled,Ppfiled,Rfiled
    bpfiled,cpfiled     
    cdinfofiled
    x0filed,y0filed,z0filed
    alphafiled
*/

  // Load matrix A
   ierr = MatLoadTransposefromFile(PETSC_COMM_WORLD,"Apfiled",MATMPIAIJ,&A);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix A loaded.\n");CHKERRQ(ierr);  
  // Load the preconditioner matrix Ac
    ierr = MatLoadTransposefromFile(PETSC_COMM_WORLD,"Acfiled",MATMPIAIJ,&Ac);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix Ac received.\n");CHKERRQ(ierr);  
  // Load the mapping matrix Prt (the transepose of Pr)
    ierr = MatLoadTransposefromFile(PETSC_COMM_WORLD,"Ppfiled",MATMPIAIJ,&Prt);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix Prt received.\n");CHKERRQ(ierr);  
  // Load the permutation matrix R 
    ierr = MatLoadTransposefromFile(PETSC_COMM_WORLD,"Rfiled",MATMPIAIJ,&R);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: matrix R received.\n");CHKERRQ(ierr);  
  // Load RHS vector b
    ierr = VecLoadfromFile(PETSC_COMM_WORLD,"bpfiled",VECMPI,&b);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: vector b received.\n");CHKERRQ(ierr);   
  // Load RHS vector c 
    ierr = VecLoadfromFile(PETSC_COMM_WORLD,"cpfiled",VECMPI,&c);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: vector c received.\n");CHKERRQ(ierr);  

  // Load cdinfo
    ierr = CdinfoLoad(PETSC_COMM_WORLD,"cdinfofiled", &cdinfo);CHKERRQ(ierr);
  // Load x0, y0, z0 as initial points for streamlines  
    ierr =  VecLoadfromFile(PETSC_COMM_WORLD,"x0filed",VECMPI,&x0);CHKERRQ(ierr);
    ierr =  VecLoadfromFile(PETSC_COMM_WORLD,"y0filed",VECMPI,&y0);CHKERRQ(ierr);
    ierr =  VecLoadfromFile(PETSC_COMM_WORLD,"z0filed",VECMPI,&z0);CHKERRQ(ierr);

  // Load betainp, this is the design parameters
//    ierr = MatLoadTransposefromFile(PETSC_COMM_WORLD,"alphafile",MATMPIAIJ,&alphamat);CHKERRQ(ierr);
  
   if(restart==0){ierr = VecLoadfromFile(PETSC_COMM_WORLD,"beinpfiled",VECMPI,&beinp);CHKERRQ(ierr);}
             else{ierr = VecLoadfromFile(PETSC_COMM_WORLD,"alphabinfiles",VECMPI,&beinp);CHKERRQ(ierr);}

   if(region==1){  ierr = VecLoadfromFile(PETSC_COMM_WORLD,"beregionfiled",VECMPI,&regionvec);CHKERRQ(ierr);}      

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//Preprocess loaded data 

  // Form x,y and adjvec by copying b       
          ierr = VecDuplicate(b,&x);CHKERRQ(ierr);
          ierr = VecDuplicate(b,&y);CHKERRQ(ierr);
          ierr = VecDuplicate(b,&beinA);CHKERRQ(ierr);   
          ierr = VecDuplicate(b,&adjvec);CHKERRQ(ierr); 
      
          ierr = VecDuplicate(c,&cp);CHKERRQ(ierr);

  // Form Ap by copying A    
          ierr = MatDuplicate(A,MAT_COPY_VALUES,&Ap);CHKERRQ(ierr);
  // Get some number 
          // m: number of x
          // n: number of alpha
          ierr = MatGetSize(Prt,&m,&n);CHKERRQ(ierr);

  // Set vector gradf 
	  ierr = VecCreate(PETSC_COMM_WORLD,&gradf);CHKERRQ(ierr);
      	  ierr = VecSetSizes(gradf,PETSC_DECIDE,m);CHKERRQ(ierr);
      	  ierr = VecSetType(gradf, VECMPI);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Completed data receiving and setting!\n\n");CHKERRQ(ierr);
  // ctx scatters parallel velocity solution to sequential one
     ierr = VecScatterCreateToAll(y,&ctx,&v);CHKERRQ(ierr);


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//Presetting
  //Ksp
  	  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp1);CHKERRQ(ierr);
  	  ierr = KSPSetTolerances(ksp1,rtol,abstol,PETSC_DEFAULT,kspmaxits);CHKERRQ(ierr);
  	  ierr = KSPSetFromOptions(ksp1);CHKERRQ(ierr);
          ierr = KSPSetInitialGuessNonzero(ksp1,PETSC_FALSE);CHKERRQ(ierr);
          
  	  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp2);CHKERRQ(ierr);
          ierr = KSPSetTolerances(ksp2,rtol,abstol,PETSC_DEFAULT,kspmaxits);CHKERRQ(ierr);
  	  ierr = KSPSetFromOptions(ksp2);CHKERRQ(ierr);
          ierr = KSPSetInitialGuessNonzero(ksp2,PETSC_FALSE);CHKERRQ(ierr);
  //Eps
          // Can't set it here!
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// Begin optimization      
for(iter=0;iter<niter;iter++){            
      // Set alphamat 

         ierr =  MatMultTranspose(Prt,beinp,beinA);CHKERRQ(ierr);
         ierr =  MatDiagonalSet(Ap,beinA,ADD_VALUES);CHKERRQ(ierr);        
  
         ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Alphavec is received and set \n");CHKERRQ(ierr);
      
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// Solve Velocity          


         ierr = PetscGetTime(&t1);CHKERRQ(ierr);
  	 ierr = KSPSetOperators(ksp1,Ap,Ac,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  	 ierr = KSPSolve(ksp1,b,x);CHKERRQ(ierr);
  	 ierr = PetscGetTime(&t2);CHKERRQ(ierr);
 	 elapsed_time = t2 - t1;
  	 ierr = KSPGetIterationNumber(ksp1,&its);CHKERRQ(ierr); 
  	 PetscPrintf(PETSC_COMM_WORLD,"PETSC: Velocity field is solved in %4.2e secs with %d iterations \n",elapsed_time,its);
  

         // Set the current solution to be the initial guess of next iteraiton
         ierr = KSPSetInitialGuessNonzero(ksp1,PETSC_TRUE);CHKERRQ(ierr);
         
 
         //scale velocity and pressure to a small number
         // PetscScalar  xmax;
         // VecMax(x,PETSC_NULL,&xmax);
         // VecScale(x,0.01/xmax);
         // VecScale(b,0.01/xmax);
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Reverse permutation
          ierr = MatMult(R,x,y);CHKERRQ(ierr);
          // Now y is the unpermuted order velocity and pressure field
          
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
//  Save x to a binary file  

          //binary file velocity and pressure
          //ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"solbin",FILE_MODE_WRITE,&binv);CHKERRQ(ierr);
          //ierr = VecView(y,  binv);CHKERRQ(ierr);
          //ierr = PetscViewerDestroy(binv);CHKERRQ(ierr);
          //ascii file velocity and pressure 
          //ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"sol",&ascv);CHKERRQ(ierr);
          //ierr = VecView(y,  ascv);CHKERRQ(ierr);
          //ierr = PetscViewerDestroy(ascv);CHKERRQ(ierr);

          ierr = VecSave(y,"solfiles");CHKERRQ(ierr);
          ierr = VecSavebin(y,"solbinfiles");CHKERRQ(ierr);

          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: v and gradf have been stored.\n\n");CHKERRQ(ierr);

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Now scatter velocity solution to each process, prepare for streamline integration 
          ierr = VecScatterBegin(y,v,INSERT_VALUES,SCATTER_FORWARD,ctx);CHKERRQ(ierr);  
          ierr = VecScatterEnd(y,v,INSERT_VALUES,SCATTER_FORWARD,ctx);CHKERRQ(ierr); 
          if(withMatlab==1){VecView(y,socketviewer);}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Solve Markov Map  

          if(operator==0){ MarkovMap2d(v, ny, nz, cdinfo,&MarkovA,&dAdy,&dAdz);}
                     else{ MarkovMap2dv(v, ny, nz, cdinfo,&MarkovA,&dAdy,&dAdz);}  
           
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Markov matrix is obtained, now saving .\n\n");CHKERRQ(ierr);       
          SparseSave(MarkovA,"A");

          
          if(withMatlab==1){SparseView(MarkovA,socketviewer);}
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
// Eigenvalue of Markov Matrix

      if(eig==1){     
          ierr = MatGetVecs(MarkovA,PETSC_NULL,&xr);CHKERRQ(ierr);
          ierr = MatGetVecs(MarkovA,PETSC_NULL,&xi);CHKERRQ(ierr);

          ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
          ierr = EPSSetOperators(eps,MarkovA,PETSC_NULL);CHKERRQ(ierr);
          ierr = EPSSetProblemType(eps,problemtype);CHKERRQ(ierr);
          ierr = EPSSetType(eps,type);
          ierr = EPSSetDimensions(eps,nev,ncv);
          ierr = EPSSetTolerances(eps,epstol,maxit);CHKERRQ(ierr);
          ierr = EPSSetWhichEigenpairs(eps,EPS_LARGEST_REAL);CHKERRQ(ierr);        

          ierr = EPSSolve(eps);CHKERRQ(ierr);
          ierr = EPSGetIterationNumber(eps, &its);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Number of iterations of the method: %d\n",its);CHKERRQ(ierr);

       	  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  	  ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Number of converged eigenpairs: %d\n\n",nconv);CHKERRQ(ierr);
   	  
   	  ierr = EPSGetDimensions(eps,&nev,PETSC_NULL);CHKERRQ(ierr);
  
          for(i=0;i<nev;i++){
    	   ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr); 
    	   ierr = PetscPrintf(PETSC_COMM_WORLD,"%f +i %f\n",kr,ki);CHKERRQ(ierr);
          }
          //get the SLEV
          ierr = EPSGetEigenpair(eps,1,&kr,&ki,xr,xi);CHKERRQ(ierr); 
       }
///////////////////////////////////////////////////////////////////////////
switch (app){
   case 1: // given c (permutated)
          // be careful, c has to be permutated
          break;
   case 2: // dfdv is calculated by given function
          // first zero c!
          ierr  = VecScale (c, 0);CHKERRQ(ierr);

          dMarkovMap2df(v, ny, nz, cdinfo,&dldv1);
          // scatter size(v) to size(sol)
          ierr  = VecGetOwnershipRange(dldv1,&Istart,&Iend);CHKERRQ(ierr);
          Lsize = Iend-Istart;
          ierr  = ISCreateStride(MPI_COMM_WORLD,Lsize,Istart,1,&is2);CHKERRQ(ierr);
          ierr  = VecScatterCreate(dldv1,is2,c,is2,&ctx1);CHKERRQ(ierr);
          ierr  = VecScatterBegin(dldv1,c,INSERT_VALUES,SCATTER_FORWARD,ctx1);CHKERRQ(ierr);
          ierr  = VecScatterEnd(dldv1,c,INSERT_VALUES,SCATTER_FORWARD,ctx1);CHKERRQ(ierr);
          // Permutation
          ierr = MatMultTranspose(R,c,cp);CHKERRQ(ierr);
          ierr = VecCopy(cp,c);CHKERRQ(ierr);

          break;

   case 3: // second largest eigenvalue
         // first zero c!
         ierr  = VecScale (c, 0);CHKERRQ(ierr);
         dMarkovMap2d(v, ny, nz, cdinfo,&dydv, &dzdv);

         // Solve dlambda/dalpha
         ierr  = MatGetVecs(MarkovA,PETSC_NULL,&Ones);CHKERRQ(ierr);
    	 ierr  = MatGetVecs(MarkovA,PETSC_NULL,&dldy);CHKERRQ(ierr);
    	 ierr  = MatGetVecs(MarkovA,PETSC_NULL,&dldz);CHKERRQ(ierr);
    	 ierr  = MatGetVecs(dydv,&dldv1,PETSC_NULL);CHKERRQ(ierr);
    	 ierr  = MatGetVecs(dzdv,&dldv2,PETSC_NULL);CHKERRQ(ierr);

    	 ierr  = VecGetOwnershipRange(Ones,&Istart,&Iend);CHKERRQ(ierr);
         Lsize = Iend-Istart;
    	 for(i=Istart;i<Iend;i++) {VecSetValue(Ones,i,1,INSERT_VALUES );}
    
    	 
    	 ierr  = VecNorm(xr,2,&xrnorm);CHKERRQ(ierr);  
    	 ierr  = VecScale(xr, 1.0/xrnorm);CHKERRQ(ierr);

    	 ierr  = MatDiagonalScale(dAdy,xr,xr);CHKERRQ(ierr);
    	 ierr  = MatDiagonalScale(dAdz,xr,xr);CHKERRQ(ierr);

    	 ierr  = MatMult(dAdy,Ones,dldy);CHKERRQ(ierr);
    	 ierr  = MatMult(dAdz,Ones,dldz);CHKERRQ(ierr);

	 ierr  = MatMultTranspose(dydv,dldy,dldv1);CHKERRQ(ierr);
    	 ierr  = MatMultTranspose(dzdv,dldz,dldv2);CHKERRQ(ierr);

    	 ierr  = VecAXPY(dldv1,1,dldv2);CHKERRQ(ierr);

    	 ierr  = VecGetSize(dldv1,&dldvlength);CHKERRQ(ierr);
    	 ierr  = VecGetOwnershipRange(c,&Istart,&Iend);CHKERRQ(ierr);
    	 Lsize = Iend-Istart;
    
          // scatter size(v) to size(sol)
          ierr  = VecGetOwnershipRange(dldv1,&Istart,&Iend);CHKERRQ(ierr);
          Lsize = Iend-Istart;
          ierr  = ISCreateStride(MPI_COMM_WORLD,Lsize,Istart,1,&is2);CHKERRQ(ierr);
          ierr  = VecScatterCreate(dldv1,is2,c,is2,&ctx1);CHKERRQ(ierr);
          ierr  = VecScatterBegin(dldv1,c,INSERT_VALUES,SCATTER_FORWARD,ctx1);CHKERRQ(ierr);
          ierr  = VecScatterEnd(dldv1,c,INSERT_VALUES,SCATTER_FORWARD,ctx1);CHKERRQ(ierr);

          ierr  = ISDestroy(is2);CHKERRQ(ierr);
          ierr  = VecScatterDestroy(ctx1);CHKERRQ(ierr);
          // Permutation
          ierr = MatMultTranspose(R,c,cp);CHKERRQ(ierr);
          ierr = VecCopy(cp,c);CHKERRQ(ierr);
          ierr = VecScale (c, -1);CHKERRQ(ierr);
          break;

}



//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Solve adjoint equation

     ierr = PetscGetTime(&t1);CHKERRQ(ierr);
     ierr = KSPSetOperators(ksp2,Ap,Ac,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
     ierr = KSPSolve(ksp2,c,adjvec);CHKERRQ(ierr);
     ierr = PetscGetTime(&t2);CHKERRQ(ierr);
     elapsed_time = t2 - t1;
     ierr = KSPGetIterationNumber(ksp2,&its);CHKERRQ(ierr);  
     PetscPrintf(PETSC_COMM_WORLD,"PETSC: Adjoint equation is solved in %4.2e secs with %d iterations \n",elapsed_time,its);


     // Recover Ap to be A
     ierr = MatDuplicate(A,MAT_COPY_VALUES,&Ap);CHKERRQ(ierr);
     // Set the current solution to be the initial guess of next iteraiton
     ierr = KSPSetInitialGuessNonzero(ksp2,PETSC_TRUE);CHKERRQ(ierr);

// Postprocess to find gradient vector and fobj  
          // fobj
  	  ierr = VecDot(c,x,&fobj);CHKERRQ(ierr);
  	  // gradf
  	  ierr = VecPointwiseMult(y, x,adjvec);CHKERRQ(ierr);
  	  ierr = MatMult(Prt,y,gradf);CHKERRQ(ierr);


          // Here we constrain the material on top and bottom
          //ierr  = VecGetOwnershipRange(gradf,&Istart,&Iend);CHKERRQ(ierr);
          //Lsize = Iend-Istart;
          //ierr  = VecGetArray(gradf,&gradfarray);CHKERRQ(ierr);
          //if(rank>1&&rank<size-2){for(i=0;i<Lsize;i++){*(gradfarray+i)=0;}}
          //ierr  = VecRestoreArray(gradf,&gradfarray);CHKERRQ(ierr);
          


          ierr = VecNorm(gradf,2,&gradfnorm);CHKERRQ(ierr);  
          ierr = VecScale (gradf, 1.0/gradfnorm);CHKERRQ(ierr);  
          ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: fobj and gradf have been found and normalized.\n");CHKERRQ(ierr);



          ierr = VecSave(gradf,"gradffiles");CHKERRQ(ierr); 
          /*
          //binary file gradf
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"gradfbin",FILE_MODE_WRITE,&binv);CHKERRQ(ierr);
          ierr = VecView(gradf,  binv);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(binv);CHKERRQ(ierr);
          //ascii file gradf 
          ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,"gradf",&ascv);CHKERRQ(ierr);
          ierr = VecView(gradf,  ascv);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(ascv);CHKERRQ(ierr);
          */


          //if(withMatlab==1){VecView(gradf,socketviewer);}

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// alpha updating
    
    //new scale scheme
    PetscInt maxgradfnum;

    if(region==1) { VecPointwiseMult(gradf,gradf,regionvec);}
    ierr = VecMax(gradf,&maxgradfnum,&gradfnorm);CHKERRQ(ierr);
    gradfnorm  =stepsize*1.0/gradfnorm;

    //stepsize =stepsize*1.1;
    ierr = VecAXPY(beinp,gradfnorm,gradf);CHKERRQ(ierr);
     
     ierr = PetscPrintf(PETSC_COMM_WORLD,"PETDC: The new step size %f .\n",gradfnorm);CHKERRQ(ierr);
  
     //ierr = VecAXPY(beinp,stepsize,gradf);CHKERRQ(ierr);
  

   ierr  = VecGetArray(beinp,&beinparray);CHKERRQ(ierr);
   ierr  = VecGetOwnershipRange(beinp,&Istart,&Iend);CHKERRQ(ierr);
   Lsize = Iend-Istart;

   //set the negative part of beinparray to be 0
   for(i=0;i<Lsize;i++) {if(*(beinparray+i)<0){*(beinparray+i)=0;} }
   // optionally make the sum to be 0 (fix material used) 
   ierr  = VecRestoreArray(beinp,&beinparray);CHKERRQ(ierr);

    if(region==1) { VecPointwiseMult(beinp,beinp,regionvec);}


   sprintf(alphafilename,"alpha%dfiles",iter);
   
   ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: File %s is saved.\n",alphafilename);CHKERRQ(ierr);

   ierr = VecSave(beinp,alphafilename);CHKERRQ(ierr);
   ierr = VecSavebin(beinp,"alphabinfiles");CHKERRQ(ierr);

if(withMatlab==1){VecView(beinp,socketviewer);}
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Destroy Everything       


  PetscPrintf(PETSC_COMM_WORLD,"PETSC:Destroy variables\n\n");
 
  //Destroy: V solve
  ierr = KSPDestroy(ksp1);CHKERRQ(ierr);
  ierr = KSPDestroy(ksp2);CHKERRQ(ierr);
  ierr = VecDestroy(b);CHKERRQ(ierr);
  ierr = VecDestroy(x);CHKERRQ(ierr);
  ierr = VecDestroy(y);CHKERRQ(ierr);
  ierr = VecDestroy(c);CHKERRQ(ierr);
  if(app==2||app==3){ierr = VecDestroy(cp);CHKERRQ(ierr);}
  ierr = VecDestroy(beinp);CHKERRQ(ierr);
  ierr = VecDestroy(beinA);CHKERRQ(ierr);
  ierr = VecDestroy(adjvec);CHKERRQ(ierr);
  ierr = MatDestroy(A);CHKERRQ(ierr);
  ierr = MatDestroy(Ac);CHKERRQ(ierr);
  ierr = MatDestroy(Ap);CHKERRQ(ierr); // ierr = MatDestroy(alphamat);CHKERRQ(ierr);
  //Destroy: Markov Matrix  
    
  ierr = MatDestroy(MarkovA);CHKERRQ(ierr);
  ierr = MatDestroy(dAdy);CHKERRQ(ierr);
  ierr = MatDestroy(dAdz);CHKERRQ(ierr);
  if(app==3)
  {
   ierr = MatDestroy(dydv);CHKERRQ(ierr);
   ierr = MatDestroy(dzdv);CHKERRQ(ierr); 
  }

   ierr = VecScatterDestroy(ctx);CHKERRQ(ierr);
   ierr = VecDestroy(v);CHKERRQ(ierr);
 //  ierr = VecDestroy(x0);CHKERRQ(ierr);
 //  ierr = VecDestroy(z0);CHKERRQ(ierr);
 //  ierr = VecDestroy(y0);CHKERRQ(ierr);
 
  //Destroy: Eigenvalue
  if(eig==1){
  ierr = EPSDestroy(eps);CHKERRQ(ierr);
  }
  //Destroy: Gradient!
  if(app==3)
  {
  ierr = VecDestroy(dldy);CHKERRQ(ierr);
  ierr = VecDestroy(dldz);CHKERRQ(ierr);
  ierr = VecDestroy(dldv2);CHKERRQ(ierr);
  }
  if(app==2||app==3) {ierr = VecDestroy(dldv1);CHKERRQ(ierr);}

  //ierr = PetscViewerDestroy(solfileviewer);CHKERRQ(ierr);
  //ierr = PetscViewerDestroy(markovfileviewer);CHKERRQ(ierr);

  //ierr = PetscFinalize();CHKERRQ(ierr); 
  ierr = SlepcFinalize();CHKERRQ(ierr);
  return 0;
}



