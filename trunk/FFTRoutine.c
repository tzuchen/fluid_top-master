#include "petscksp.h"
#include "petscda.h"
#include </home/tzuchen/fftw/include/fftw3.h>
//#include <fftw3.h> 
#include <math.h>
#include <stdlib.h>
#include "FFTRoutine.h"
#include "LargeVecFunction.h"
//////////////////////////////////////////////////////////////////////
int FFT2D(Vec x, Vec y, PetscInt m, PetscInt n, PetscInt cIstart, PetscInt cIend, PetscInt iter, PetscTruth Smoothingflag, PetscScalar diff, char fnameadd[])
{
  
   PetscErrorCode   ierr;
   PetscInt         i,wnmax,rank,ncol;
   PetscScalar      *Erecord;
   FILE             *fid;  
   char             wnfname[50];
   PetscLogDouble v1,v2,elapsed_time;

   ierr = PetscPrintf(PETSC_COMM_WORLD,"\nPETSC: Doing 2D FFT Now! \n");
   MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
   wnmax    = floor(n*0.72);
   Erecord  = malloc(wnmax*sizeof(PetscScalar));
   ncol     = cIend - cIstart;

   
 

   ierr = PetscGetTime(&v1);CHKERRQ(ierr);
       SkewConvert(x,y,(PetscInt)(n*0.5),cIstart,cIend);
   ierr = PetscGetTime(&v2);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: SkewConverting takes %f sec\n",v2-v1);

   ierr = PetscGetTime(&v1);CHKERRQ(ierr);
       FFTColumn(x,y,ncol,n);


   ierr = PetscGetTime(&v2);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Column FFT takes %f sec\n",v2-v1);

   ierr = PetscGetTime(&v1);CHKERRQ(ierr);
     //  FFTRow(x,y, (PetscInt)(n*0.5), n,cIstart,cIend,&Erecord, wnmax,diff);
       FFTRow3(x,y, (PetscInt)(n*0.5), n,cIstart,cIend,&Erecord, wnmax,diff);
   ierr = PetscGetTime(&v2);CHKERRQ(ierr);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Row FFT takes %f sec\n",v2-v1);

       IFFTColumn(x,y,ncol,n);
   ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: Column IFFT done\n");

   if(Smoothingflag==PETSC_FALSE) {sprintf(wnfname, "wnlog_%d_%d%s",n,iter,fnameadd);}
                  else{sprintf(wnfname, "wnlog_%d_%d%ss",n,iter,fnameadd);}

   if(rank==0){   
     fid = fopen(wnfname,"w");
     for(i=0;i<wnmax;i++){fprintf(fid,"%d   %f \n", i,*(Erecord+i));}
     fclose(fid);
   }
  
   free(Erecord);

 return 0;
}
//////////////////////////////////////////////////////////////////////
/* Vec x and y have the same parallel layout 
   Vec x is the data, this function writes vec y as 
   y(i) = x(n-i);
 
*/

int SkewConvert(Vec x, Vec y, PetscInt nrow,PetscInt cIstart, PetscInt cIend)
{
   
    PetscErrorCode  ierr;
    PetscInt        size,rank;
    PetscInt        colcount,i,j,ncol;  
    VecScatter      ctx;
    PetscInt        pmax,puse,k; 
    PetscScalar     *xpt,*ypt,*ylpt; 
    PetscInt        *indI,*indJ;
    Vec             yl;

    PetscInt        nvec=2;
    Vec             xvec[nvec];
    IS              ISfrom[nvec],ISto[nvec];
    VecScatter      ctxt[nvec];             

 
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    pmax = 1e6;
    ncol = (PetscInt)(pmax*1.0/nrow-1);


    indI = malloc(pmax*sizeof(PetscInt));
    indJ = malloc(pmax*sizeof(PetscInt));
   

    VecGetArray(y,&ypt);
    LargeVecCreate(&x,nvec,xvec);
    colcount  = cIstart;


 while(colcount<cIend){

     if(colcount+ncol>=cIend){ncol = cIend-colcount;}
     puse = ncol*nrow;
     VecCreateSeqWithArray(MPI_COMM_SELF,puse,ypt,&yl);

     k = 0;  
     for(i=colcount;i<colcount+ncol;i++){
         for(j=0;j<nrow;j++){
            *(indI+k) = 2*nrow-i-1;
            *(indJ+k) = nrow-j-1; 

            k++;

         }
     }  
 
    ISCreateGeneralWithIJ(MPI_COMM_SELF,x,xvec,nvec,nrow,puse, indI, indJ,ISfrom, ISto);
    LargeVecScatterCreate(xvec,ISfrom,yl,ISto ,ctxt,nvec); 
    ISArrayDestroy(ISfrom,nvec);
    ISArrayDestroy(ISto,nvec);
    LargeVecScatterBeginEnd(xvec,yl,INSERT_VALUES,SCATTER_FORWARD,ctxt,nvec);
    VecScatterArrayDestroy(ctxt,nvec);

    VecDestroy(yl);
    ypt = ypt+puse;
    colcount = colcount+ncol;
  

  }


free(indI);
free(indJ);


return 0;
}


///////////////////////////////////////////////////////////////////////
/*
This function do FFT on columns of the matrix [x;y]
it stores the real part in x, and the image pare in y; 
*/     

int FFTColumn(Vec x, Vec y, PetscInt m, PetscInt n){
  
   PetscScalar    *xpt,*ypt;
   PetscInt       i,k,n2;
   fftw_complex   *in, *out;
   fftw_plan      p;

  
   n2  = (PetscInt)(n*0.5);
   in  = fftw_malloc(n*sizeof(fftw_complex));
   out = fftw_malloc(n*sizeof(fftw_complex));


   VecGetArray(x,&xpt);
   VecGetArray(y,&ypt);

   for(k=0;k<m;k++) {
        
       for(i=0;i<n2;i++){  in[i][0] = *(xpt+i); in[i][1]= 0; }  
       for(i=0;i<n2;i++){  in[i+n2][0] = *(ypt+i); in[i+n2][1]= 0; }  
  
         p = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);         
         fftw_execute(p);
         fftw_destroy_plan(p);
        // fprintf(stdout, "Doing FFT m=%d %d \n",m,k);
 
       for(i=0;i<n2;i++){ *(xpt+i) = out[i][0]; }  
       for(i=0;i<n2;i++){ *(ypt+i) = out[i][1]; } 
 
         xpt = xpt+ n2;
         ypt = ypt+ n2;
    }
 


      fftw_free(in); fftw_free(out);

      VecRestoreArray(x,&xpt);
      VecRestoreArray(y,&ypt);

      return 0;
   }


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


int FFTRow(Vec x, Vec y, PetscInt grow, PetscInt gcol, PetscInt cIstart, PetscInt cIend, PetscScalar **Erecordpt, PetscInt wnmax,PetscScalar diff)
{

   PetscErrorCode  ierr;
   PetscInt        rank,size;
   PetscInt        i,j,k;
   PetscInt        *Ecount,wn;
   PetscInt        ncol,nrow;
   PetscScalar     *Erecord,zf;
   PetscScalar     *xpt,*ypt;
   fftw_complex    *in, *out;
   fftw_plan       p;
   MPI_Status      status;
   PetscScalar     *temprowRe,*temprowIm,*Rept,*Impt;
   PetscInt        rowcount;  


    ncol = cIend-cIstart;
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

   PetscInt        ncolarray[size],disp[size];
   
   for(i=0;i<size;i++){
   MPI_Gather(&ncol,1, MPI_INT, ncolarray, 1,  MPI_INT, i,PETSC_COMM_WORLD);
   }
   disp[0] = 0;//ncolarray[0];
   for(i=1;i<size;i++){
   disp[i] = disp[i-1]+ncolarray[i-1];
   }


   VecGetArray(x,&xpt);
   VecGetArray(y,&ypt);
   Erecord  = *Erecordpt;
   Ecount   = malloc(wnmax*sizeof(PetscInt)); 

   for(i=0;i<wnmax;i++){*(Ecount+i) = 0;*(Erecord+i) = 0;}
   in  = fftw_malloc(gcol*sizeof(fftw_complex));
   out = fftw_malloc(gcol*sizeof(fftw_complex));



   temprowRe =  malloc(ncol*sizeof(PetscScalar));
   temprowIm =  malloc(ncol*sizeof(PetscScalar));
   Rept      = malloc(gcol*sizeof(PetscScalar));
   Impt      = malloc(gcol*sizeof(PetscScalar));
   
  
   
   rowcount = 0 ;
while(rowcount<grow){

    
     for(j=0;j<size;j++){
         for(i=0;i<ncol;i++){ 
             if(rowcount+j<grow){
                *(temprowRe+i) = *(xpt+grow*i+rowcount+j);
                *(temprowIm+i) = *(ypt+grow*i+rowcount+j); 
             }
         }
         MPI_Gatherv(temprowRe, ncol, MPI_DOUBLE, Rept, ncolarray, disp, MPI_DOUBLE, j, PETSC_COMM_WORLD);
         MPI_Gatherv(temprowIm, ncol, MPI_DOUBLE, Impt, ncolarray, disp, MPI_DOUBLE, j, PETSC_COMM_WORLD);
   }   
   if(rowcount+rank<grow){
      for(i=0;i<gcol;i++){  in[i][0] = *(Rept+i); in[i][1]= *(Impt+i); } 

      p = fftw_plan_dft_1d(gcol, in, out, FFTW_FORWARD, FFTW_ESTIMATE);         
      fftw_execute(p);
      fftw_destroy_plan(p); 
      for(i=0;i<gcol;i++){

            wn = WaveNumber(i,rowcount+rank,gcol);
            zf=1;
            *(Erecord+wn) = *(Erecord+wn) + zf*sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]);
            *(Ecount+wn)  = *(Ecount+wn)+1; 
            // Do diffusion
            out[i][0] = out[i][0]*exp(-4*M_PI*M_PI*wn*wn*diff);
            out[i][1] = out[i][1]*exp(-4*M_PI*M_PI*wn*wn*diff);
      }
      // Do IFFT
      p = fftw_plan_dft_1d(gcol, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);         
      fftw_execute(p);
      fftw_destroy_plan(p); 
      for(i=0;i<gcol;i++){*(Rept+i) = in[i][0]/gcol ; *(Impt+i) = in[i][1]/gcol; }

   }  
     
    //Scatter Backward
    for(j=0;j<size;j++){ 
     MPI_Scatterv ( Rept, ncolarray,disp, MPI_DOUBLE, temprowRe, ncol, MPI_DOUBLE,j, PETSC_COMM_WORLD); 
     MPI_Scatterv ( Impt, ncolarray,disp, MPI_DOUBLE, temprowIm, ncol, MPI_DOUBLE,j, PETSC_COMM_WORLD); 
         for(i=0;i<ncol;i++){ 
             if(rowcount+j<grow){
                 *(xpt+grow*i+rowcount+j) = *(temprowRe+i);
                 *(ypt+grow*i+rowcount+j) = *(temprowIm+i); 
             }
         }
     }



   rowcount = rowcount+size;
}  

  
  
         if(rank==0){
       PetscScalar    *Erecord2; 
       PetscInt       *Ecount2;  
       //MPI_Status     status;
       
       Erecord2  = malloc(wnmax*sizeof(PetscScalar));
       Ecount2   = malloc(wnmax*sizeof(PetscInt));
      
         for(i=1;i<size;i++){
            MPI_Recv(Erecord2,wnmax,MPI_DOUBLE,i,i,PETSC_COMM_WORLD,&status);
            MPI_Recv(Ecount2,wnmax,MPI_INT,i,i*10,PETSC_COMM_WORLD,&status);
            for(j=0;j<wnmax;j++){
                  *(Erecord+j) = *(Erecord+j) + *(Erecord2+j);
                  *(Ecount+j) = *(Ecount+j) + *(Ecount2+j);
            }  
         } 
           free(Erecord2);
           free(Ecount2);       
       }else{
        MPI_Send(Erecord,wnmax,MPI_DOUBLE,0,rank,PETSC_COMM_WORLD);
        MPI_Send(Ecount,wnmax,MPI_INT,0,rank*10,PETSC_COMM_WORLD);
      } 

      //normalize
      if(rank==0){ for(i=0;i<wnmax;i++){*(Erecord+i) = *(Erecord+i)/gcol*1.0/gcol;}}

  





   free(Rept);
   free(Impt);
   free(temprowRe);
   free(temprowIm);
   free(Ecount);

 return 0;

}


//////////////////////////////////////////////////////////////////////////////////////

int FFTRow2(Vec x, Vec y, PetscInt grow, PetscInt gcol, PetscInt cIstart, PetscInt cIend, PetscScalar **Erecordpt, PetscInt wnmax,PetscScalar diff)
{

   PetscErrorCode  ierr;
   PetscInt        rank,size;
   PetscInt        i,j,k;
   PetscInt        *Ecount,wn;
   PetscInt        ncol,nrow;
   PetscScalar     *Erecord,zf;
   PetscScalar     *xpt,*ypt;
   fftw_complex    *in, *out;
   fftw_plan       p;
   MPI_Status      status;
   PetscScalar     *temprowRe,*temprowIm,*Rept,*Impt;
   PetscInt        rowcount;  
   PetscInt        pmax = 1e7;           
   PetscLogDouble  t1,t2;
   MPI_Request     request;



    ncol = cIend-cIstart;
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

   PetscInt        ncolarray[size],disp[size];
   
   for(i=0;i<size;i++){
   MPI_Gather(&ncol,1, MPI_INT, ncolarray, 1,  MPI_INT, i,PETSC_COMM_WORLD);
   }



   VecGetArray(x,&xpt);
   VecGetArray(y,&ypt);
   Erecord  = *Erecordpt;
   Ecount   = malloc(wnmax*sizeof(PetscInt)); 

   for(i=0;i<wnmax;i++){*(Ecount+i) = 0;*(Erecord+i) = 0;}
   in  = fftw_malloc(gcol*sizeof(fftw_complex));
   out = fftw_malloc(gcol*sizeof(fftw_complex));

   nrow = (PetscInt)floor(pmax*1.0/gcol);
   //nrow = 4;
   if(nrow*size>grow){nrow =(PetscInt)ceil(grow*1.0/size);}

   for(i=0;i<size;i++){ncolarray[i]=ncolarray[i]*nrow;}
   disp[0] = 0;//ncolarray[0];
   for(i=1;i<size;i++){
   disp[i] = disp[i-1]+ncolarray[i-1];
   }

   temprowRe =  malloc(nrow*ncol*sizeof(PetscScalar));
   temprowIm =  malloc(nrow*ncol*sizeof(PetscScalar));
   Rept      = malloc(nrow*gcol*sizeof(PetscScalar));
   Impt      = malloc(nrow*gcol*sizeof(PetscScalar));
   
   PetscInt  st,rf;

 
   fprintf(stdout,"Begin While loop %d %d %d\n",size,nrow, grow);
   rowcount = 0 ;
while(rowcount<grow){

  ierr = PetscGetTime(&t1);CHKERRQ(ierr);
    for(j=0;j<size;j++){

       st = (rank+j)%size;
       rf = (rank-j)%size;
       if(rf<0){rf=rf+size;}        


       for(k=0;k<nrow;k++){  
          for(i=0;i<ncol;i++){ 
              if(rowcount+k*size+j<grow){
                 *(temprowRe+k*ncol+i) = *(xpt+grow*i+rowcount+k*size+st);
                 *(temprowIm+k*ncol+i) = *(ypt+grow*i+rowcount+k*size+st); 
              }
          }
       }



         MPI_Sendrecv(temprowRe,ncol*nrow,MPI_DOUBLE,st,(rank+1)*1000+(st+1)*10+1,
                      Rept+disp[rf],ncolarray[rf],MPI_DOUBLE,rf,(rf+1)*1000+(rank+1)*10+1,PETSC_COMM_WORLD,&status);    
         MPI_Sendrecv(temprowIm,ncol*nrow,MPI_DOUBLE,st,(rank+1)*1000+(st+1)*10+2,
                      Impt+disp[rf],ncolarray[rf],MPI_DOUBLE,rf,(rf+1)*1000+(rank+1)*10+2,PETSC_COMM_WORLD,&status);
         
    }   


  ierr = PetscGetTime(&t2);CHKERRQ(ierr);


   for(k=0;k<nrow;k++){
      if(rowcount+k*size+rank<grow){
         for(i=0;i<gcol;i++){  in[i][0] = *(Rept+k*gcol+i); in[i][1]= *(Impt+k*gcol+i); } 
 
         p = fftw_plan_dft_1d(gcol, in, out, FFTW_FORWARD, FFTW_ESTIMATE);         
         fftw_execute(p);
         fftw_destroy_plan(p); 
         for(i=0;i<gcol;i++){

            wn = WaveNumber(i,rowcount+rank,gcol);
            zf=1;
            *(Erecord+wn) = *(Erecord+wn) + zf*sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]);
            *(Ecount+wn)  = *(Ecount+wn)+1; 
            // Do diffusion
            out[i][0] = out[i][0]*exp(-4*M_PI*M_PI*wn*wn*diff);
            out[i][1] = out[i][1]*exp(-4*M_PI*M_PI*wn*wn*diff);
         }
      // Do IFFT
      p = fftw_plan_dft_1d(gcol, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);         
      fftw_execute(p);
      fftw_destroy_plan(p); 
      for(i=0;i<gcol;i++){*(Rept+k*ncol+i) = in[i][0]/gcol ; *(Impt+k*ncol+i) = in[i][1]/gcol; }
     }  
   }
     
    //Scatter Backward

      for(j=0;j<size;j++){

       st = (rank+j)%size;
       rf = (rank-j)%size;
       if(rf<0){rf=rf+size;}   

         MPI_Sendrecv(Rept+disp[st],ncolarray[st],MPI_DOUBLE,st,(rank+1)*1000+(st+1)*10+3,
                      temprowRe,ncol*nrow,MPI_DOUBLE,rf,(rf+1)*1000+(rank+1)*10+3,PETSC_COMM_WORLD,&status);    
         MPI_Sendrecv(Impt+disp[st],ncolarray[st],MPI_DOUBLE,st,(rank+1)*1000+(st+1)*10+4,
                      temprowIm,ncol*nrow,MPI_DOUBLE,rf,(rf+1)*1000+(rank+1)*10+4,PETSC_COMM_WORLD,&status);

 
       for(k=0;k<nrow;k++){
         for(i=0;i<ncol;i++){ 
             if(rowcount+k*size+j<grow){
                 *(xpt+grow*i+rowcount+k*size+rf) = *(temprowRe+k*ncol+i);
                 *(ypt+grow*i+rowcount+k*size+rf) = *(temprowIm+k*ncol+i); 
             }
         }
       }
      }


   rowcount = rowcount+size*nrow;
}  

  
  
         if(rank==0){
       PetscScalar    *Erecord2; 
       PetscInt       *Ecount2;  
       //MPI_Status     status;
       
       Erecord2  = malloc(wnmax*sizeof(PetscScalar));
       Ecount2   = malloc(wnmax*sizeof(PetscInt));
      
         for(i=1;i<size;i++){
            MPI_Recv(Erecord2,wnmax,MPI_DOUBLE,i,i,PETSC_COMM_WORLD,&status);
            MPI_Recv(Ecount2,wnmax,MPI_INT,i,i*10,PETSC_COMM_WORLD,&status);
            for(j=0;j<wnmax;j++){
                  *(Erecord+j) = *(Erecord+j) + *(Erecord2+j);
                  *(Ecount+j) = *(Ecount+j) + *(Ecount2+j);
            }  
         } 
           free(Erecord2);
           free(Ecount2);       
       }else{
        MPI_Send(Erecord,wnmax,MPI_DOUBLE,0,rank,PETSC_COMM_WORLD);
        MPI_Send(Ecount,wnmax,MPI_INT,0,rank*10,PETSC_COMM_WORLD);
      } 

      //normalize
      if(rank==0){ for(i=0;i<wnmax;i++){*(Erecord+i) = *(Erecord+i)/gcol*1.0/gcol;}}

  





   free(Rept);
   free(Impt);
   free(temprowRe);
   free(temprowIm);
   free(Ecount);

 return 0;

}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int FFTRow3(Vec x, Vec y, PetscInt grow, PetscInt gcol, PetscInt cIstart, PetscInt cIend, PetscScalar **Erecordpt, PetscInt wnmax,PetscScalar diff)
{

   PetscErrorCode  ierr;
   PetscInt        rank,size;
   PetscInt        i,j,k;
   PetscInt        *Ecount,wn;
   PetscInt        ncol,nrow;
   PetscScalar     *Erecord,zf;
   PetscScalar     *xpt,*ypt;
   fftw_complex    *in, *out;
   fftw_plan       p;
   MPI_Status      status;
   PetscScalar     *temprowRe,*temprowIm,*Rept,*Impt;
   PetscInt        rowcount;  
   PetscLogDouble  t1,t2;
   MPI_Request     request;


    ncol = cIend-cIstart;
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

   PetscInt        ncolarray[size],disp[size];
   
   for(i=0;i<size;i++){
   MPI_Gather(&ncol,1, MPI_INT, ncolarray, 1,  MPI_INT, i,PETSC_COMM_WORLD);
   }
   disp[0] = 0;//ncolarray[0];
   for(i=1;i<size;i++){
   disp[i] = disp[i-1]+ncolarray[i-1];
   }


   VecGetArray(x,&xpt);
   VecGetArray(y,&ypt);
   Erecord  = *Erecordpt;
   Ecount   = malloc(wnmax*sizeof(PetscInt)); 

   for(i=0;i<wnmax;i++){*(Ecount+i) = 0;*(Erecord+i) = 0;}
   in  = fftw_malloc(gcol*sizeof(fftw_complex));
   out = fftw_malloc(gcol*sizeof(fftw_complex));



   temprowRe =  malloc(ncol*sizeof(PetscScalar));
   temprowIm =  malloc(ncol*sizeof(PetscScalar));
   Rept      = malloc(gcol*sizeof(PetscScalar));
   Impt      = malloc(gcol*sizeof(PetscScalar));
   
  
   
   rowcount = 0 ;
while(rowcount<grow){

    ierr = PetscGetTime(&t2);CHKERRQ(ierr);
     for(j=0;j<size;j++){
         for(i=0;i<ncol;i++){ 
             if(rowcount+j<grow){
                *(temprowRe+i) = *(xpt+grow*i+rowcount+j);
                *(temprowIm+i) = *(ypt+grow*i+rowcount+j); 
             }
         }
      MPI_Send(temprowRe,ncol,MPI_DOUBLE,j,rank*1000+j*10+1,PETSC_COMM_WORLD);    
      MPI_Send(temprowIm,ncol,MPI_DOUBLE,j,rank*1000+j*10+2,PETSC_COMM_WORLD);

     }   
     ierr = PetscGetTime(&t2);CHKERRQ(ierr);
//PetscPrintf(PETSC_COMM_WORLD,"send time %f \n",t2-t1);
     ierr = PetscGetTime(&t1);CHKERRQ(ierr);
      for(j=0;j<size;j++){
      MPI_Recv(Rept+disp[j],ncolarray[j],MPI_DOUBLE,j,j*1000+rank*10+1,PETSC_COMM_WORLD,&status);
      MPI_Recv(Impt+disp[j],ncolarray[j],MPI_DOUBLE,j,j*1000+rank*10+2,PETSC_COMM_WORLD,&status);
      }
     ierr = PetscGetTime(&t2);CHKERRQ(ierr);
// PetscPrintf(PETSC_COMM_WORLD,"receive time %f \n",t2-t1);   

   if(rowcount+rank<grow){
      for(i=0;i<gcol;i++){  in[i][0] = *(Rept+i); in[i][1]= *(Impt+i); } 

      p = fftw_plan_dft_1d(gcol, in, out, FFTW_FORWARD, FFTW_ESTIMATE);         
      fftw_execute(p);
      fftw_destroy_plan(p); 
      for(i=0;i<gcol;i++){

            wn = WaveNumber(i,rowcount+rank,gcol);
            if(rowcount+rank==0){zf=0.5;}else{zf=1;};
            *(Erecord+wn) = *(Erecord+wn) + zf*sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]);
            *(Ecount+wn)  = *(Ecount+wn)+1; 
            // Do diffusion
            out[i][0] = out[i][0]*exp(-4*M_PI*M_PI*wn*wn*diff);
            out[i][1] = out[i][1]*exp(-4*M_PI*M_PI*wn*wn*diff);
      }
      // Do IFFT
      p = fftw_plan_dft_1d(gcol, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);         
      fftw_execute(p);
      fftw_destroy_plan(p); 
      for(i=0;i<gcol;i++){*(Rept+i) = in[i][0]/gcol ; *(Impt+i) = in[i][1]/gcol; }

   }  
     
    //Scatter Backward
      for(j=0;j<size;j++){ 
      MPI_Send(Rept+disp[j],ncolarray[j],MPI_DOUBLE,j,rank*1000+j*10+3,PETSC_COMM_WORLD);
      MPI_Send(Impt+disp[j],ncolarray[j],MPI_DOUBLE,j,rank*1000+j*10+4,PETSC_COMM_WORLD);
      }

    for(j=0;j<size;j++){ 
           MPI_Recv(temprowRe,ncol,MPI_DOUBLE,j,j*1000+rank*10+3,PETSC_COMM_WORLD,&status);
           MPI_Recv(temprowIm,ncol,MPI_DOUBLE,j,j*1000+rank*10+4,PETSC_COMM_WORLD,&status);
 
         for(i=0;i<ncol;i++){ 
             if(rowcount+j<grow){
                 *(xpt+grow*i+rowcount+j) = *(temprowRe+i);
                 *(ypt+grow*i+rowcount+j) = *(temprowIm+i); 
             }
         }
     }



   rowcount = rowcount+size;
}  

  
  
         if(rank==0){
       PetscScalar    *Erecord2; 
       PetscInt       *Ecount2;  
       //MPI_Status     status;
       
       Erecord2  = malloc(wnmax*sizeof(PetscScalar));
       Ecount2   = malloc(wnmax*sizeof(PetscInt));
      
         for(i=1;i<size;i++){
            MPI_Recv(Erecord2,wnmax,MPI_DOUBLE,i,i,PETSC_COMM_WORLD,&status);
            MPI_Recv(Ecount2,wnmax,MPI_INT,i,i*10,PETSC_COMM_WORLD,&status);
            for(j=0;j<wnmax;j++){
                  *(Erecord+j) = *(Erecord+j) + *(Erecord2+j);
                  *(Ecount+j) = *(Ecount+j) + *(Ecount2+j);
            }  
         } 
           free(Erecord2);
           free(Ecount2);       
       }else{
        MPI_Send(Erecord,wnmax,MPI_DOUBLE,0,rank,PETSC_COMM_WORLD);
        MPI_Send(Ecount,wnmax,MPI_INT,0,rank*10,PETSC_COMM_WORLD);
      } 

      //normalize
      if(rank==0){ for(i=0;i<wnmax;i++){*(Erecord+i) = *(Erecord+i)/gcol*1.0/gcol;}}

  





   free(Rept);
   free(Impt);
   free(temprowRe);
   free(temprowIm);
   free(Ecount);

 return 0;

}


//////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
int FFTRow4(Vec x, Vec y, PetscInt grow, PetscInt gcol, PetscInt cIstart, PetscInt cIend, PetscScalar **Erecordpt, PetscInt wnmax,PetscScalar diff)
{

   PetscErrorCode  ierr;
   PetscInt        rank,size;
   PetscInt        i,j,k;
   PetscInt        *Ecount,wn;
   PetscInt        ncol,nrow;
   PetscScalar     *Erecord,zf;
   PetscScalar     *xpt,*ypt;
   fftw_complex    *in, *out;
   fftw_plan       p;
   MPI_Status      status;
   PetscScalar     *temprowRe,*temprowIm,*Rept,*Impt;
   PetscInt        rowcount;  
   PetscLogDouble  t1,t2;
   MPI_Request     request;


    ncol = cIend-cIstart;
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

   PetscInt        ncolarray[size],disp[size];
   
   for(i=0;i<size;i++){
   MPI_Gather(&ncol,1, MPI_INT, ncolarray, 1,  MPI_INT, i,PETSC_COMM_WORLD);
   }
   disp[0] = 0;//ncolarray[0];
   for(i=1;i<size;i++){
   disp[i] = disp[i-1]+ncolarray[i-1];
   }


   VecGetArray(x,&xpt);
   VecGetArray(y,&ypt);
   Erecord  = *Erecordpt;
   Ecount   = malloc(wnmax*sizeof(PetscInt)); 

   for(i=0;i<wnmax;i++){*(Ecount+i) = 0;*(Erecord+i) = 0;}
   in  = fftw_malloc(gcol*sizeof(fftw_complex));
   out = fftw_malloc(gcol*sizeof(fftw_complex));



   temprowRe =  malloc(ncol*sizeof(PetscScalar));
   temprowIm =  malloc(ncol*sizeof(PetscScalar));
   Rept      = malloc(gcol*sizeof(PetscScalar));
   Impt      = malloc(gcol*sizeof(PetscScalar));
   
   PetscInt  st,rf;
   
   
   rowcount = 0 ;
while(rowcount<grow){

    ierr = PetscGetTime(&t1);CHKERRQ(ierr);
     for(j=0;j<size;j++){
         for(i=0;i<ncol;i++){ 
             if(rowcount+j<grow){
                *(temprowRe+i) = *(xpt+grow*i+rowcount+j);
                *(temprowIm+i) = *(ypt+grow*i+rowcount+j); 
             }
         }

   //   MPI_Sendrecv(temprowRe,ncol,MPI_DOUBLE,j,rank*1000+j*10+1, 
   //                Rept+disp[j],ncolarray[j],MPI_DOUBLE,j,j*1000+rank*10+1,PETSC_COMM_WORLD,&status);
   //   MPI_Sendrecv(temprowIm,ncol,MPI_DOUBLE,j,rank*1000+j*10+2, 
   //                Impt+disp[j],ncolarray[j],MPI_DOUBLE,j,j*1000+rank*10+2,PETSC_COMM_WORLD,&status); 

       st = (rank+j)%size;
       rf = (rank-j)%size;
       if(rf<0){rf=rf+size;}
  
    // PetscPrintf(PETSC_COMM_WORLD,"st= %d rf = %d \n",st,rf);

    //  MPI_Sendrecv(temprowRe,ncol,MPI_DOUBLE,st,rank*1000+st*10+1, 
    //               Rept+disp[rf],ncolarray[rf],MPI_DOUBLE,rf,rf*1000+rank*10+1,PETSC_COMM_WORLD,&status);
  

  
     }   
for(j=0;j<size;j++){
       st = (rank+j)%size;
       rf = (rank-j)%size;
       if(rf<0){rf=rf+size;}
   // MPI_Sendrecv(temprowIm,ncol,MPI_DOUBLE,st,rank*1000+st*10+2, 
   //                Impt+disp[rf],ncolarray[rf],MPI_DOUBLE,rf,rf*1000+rank*10+2,PETSC_COMM_WORLD,&status); 
   }


     ierr = PetscGetTime(&t2);CHKERRQ(ierr);
   //  PetscPrintf(PETSC_COMM_WORLD,"send time %f \n",t2-t1);




   if(rowcount+rank<grow){
      for(i=0;i<gcol;i++){  in[i][0] = *(Rept+i); in[i][1]= *(Impt+i); } 

      p = fftw_plan_dft_1d(gcol, in, out, FFTW_FORWARD, FFTW_ESTIMATE);         
      fftw_execute(p);
      fftw_destroy_plan(p); 
      for(i=0;i<gcol;i++){

            wn = WaveNumber(i,rowcount+rank,gcol);
            zf=1;
            *(Erecord+wn) = *(Erecord+wn) + zf*sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]);
            *(Ecount+wn)  = *(Ecount+wn)+1; 
            // Do diffusion
            out[i][0] = out[i][0]*exp(-4*M_PI*M_PI*wn*wn*diff);
            out[i][1] = out[i][1]*exp(-4*M_PI*M_PI*wn*wn*diff);
      }
      // Do IFFT
      p = fftw_plan_dft_1d(gcol, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);         
      fftw_execute(p);
      fftw_destroy_plan(p); 
      for(i=0;i<gcol;i++){*(Rept+i) = in[i][0]/gcol ; *(Impt+i) = in[i][1]/gcol; }

   }  
     
    //Scatter Backward

    for(j=0;j<size;j++){ 

      MPI_Sendrecv(Rept+disp[j],ncolarray[j],MPI_DOUBLE,j,rank*1000+j*10+3,
                   temprowRe,ncol,MPI_DOUBLE,j,j*1000+rank*10+3,PETSC_COMM_WORLD,&status);
      MPI_Sendrecv(Impt+disp[j],ncolarray[j],MPI_DOUBLE,j,rank*1000+j*10+4,
                   temprowIm,ncol,MPI_DOUBLE,j,j*1000+rank*10+4,PETSC_COMM_WORLD,&status);

         //  MPI_Recv(temprowRe,ncol,MPI_DOUBLE,j,j*1000+rank*10+3,PETSC_COMM_WORLD,&status);
         //  MPI_Recv(temprowIm,ncol,MPI_DOUBLE,j,j*1000+rank*10+4,PETSC_COMM_WORLD,&status);
 
         for(i=0;i<ncol;i++){ 
             if(rowcount+j<grow){
                 *(xpt+grow*i+rowcount+j) = *(temprowRe+i);
                 *(ypt+grow*i+rowcount+j) = *(temprowIm+i); 
             }
         }
     }



   rowcount = rowcount+size;
}  

  
  
         if(rank==0){
       PetscScalar    *Erecord2; 
       PetscInt       *Ecount2;  
       //MPI_Status     status;
       
       Erecord2  = malloc(wnmax*sizeof(PetscScalar));
       Ecount2   = malloc(wnmax*sizeof(PetscInt));
      
         for(i=1;i<size;i++){
            MPI_Recv(Erecord2,wnmax,MPI_DOUBLE,i,i,PETSC_COMM_WORLD,&status);
            MPI_Recv(Ecount2,wnmax,MPI_INT,i,i*10,PETSC_COMM_WORLD,&status);
            for(j=0;j<wnmax;j++){
                  *(Erecord+j) = *(Erecord+j) + *(Erecord2+j);
                  *(Ecount+j) = *(Ecount+j) + *(Ecount2+j);
            }  
         } 
           free(Erecord2);
           free(Ecount2);       
       }else{
        MPI_Send(Erecord,wnmax,MPI_DOUBLE,0,rank,PETSC_COMM_WORLD);
        MPI_Send(Ecount,wnmax,MPI_INT,0,rank*10,PETSC_COMM_WORLD);
      } 

      //normalize
      if(rank==0){ for(i=0;i<wnmax;i++){*(Erecord+i) = *(Erecord+i)/gcol*1.0/gcol;}}

  





   free(Rept);
   free(Impt);
   free(temprowRe);
   free(temprowIm);
   free(Ecount);

 return 0;

}


//////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
/*
This function do IFFT on columns of the matrix [x;y]
 
*/     

int IFFTColumn(Vec x, Vec y, PetscInt m, PetscInt n){
  
   PetscScalar    *xpt,*ypt;
   PetscInt       i,k,n2;
   fftw_complex   *in, *out;
   fftw_plan      p;

  
   n2  = (PetscInt)(n*0.5);
   in  = fftw_malloc(n*sizeof(fftw_complex));
   out = fftw_malloc(n*sizeof(fftw_complex));


   VecGetArray(x,&xpt);
   VecGetArray(y,&ypt);

   for(k=0;k<m;k++) {
        
       for(i=0;i<n2;i++){  in[i][0] = *(xpt+i); in[i][1]= *(ypt+i); }  
       for(i=1;i<n2;i++){  in[n2+i][0] = in[n2-i][0] ; in[n2+i][1]= -in[n2-i][1]; } 

       p = fftw_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);         
       fftw_execute(p);
       fftw_destroy_plan(p);
 
       for(i=0;i<n2;i++){ *(xpt+i) = out[i][0]/n; }  
       for(i=0;i<n2;i++){ *(ypt+i) = out[i][0]/n; } 
 
         xpt = xpt+ n2;
         ypt = ypt+ n2;
    }
 
      fftw_free(in); fftw_free(out);

      VecRestoreArray(x,&xpt);
      VecRestoreArray(y,&ypt);

      return 0;
   }


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/*
 Given i,j, calculate wavenumber
*/
int WaveNumber(i,j,n){
 
   PetscInt     im,jm,wn;
   PetscScalar  ir,jr;

    if (i<n-i){im = i;}else{im=n-i;}   
    if (j<n-j){jm = j;}else{jm=n-j;}
    
  //  wn = floor(sqrt(im*im + jm*jm));

       ir = (PetscScalar)im;
       jr = (PetscScalar)jm;
  if(im==jm){wn =(PetscInt)floor(ir*sqrt(2));}else{
 
    if(jr>ir){
      wn = (PetscInt)floor(jr/cos(atan(ir/jr)));   
    }else{
      wn = (PetscInt)floor(ir/cos(atan(jr/ir))); 
    }
  }

    return wn;
}



 
