#include "petscksp.h"
#include "petscda.h"
#include <stdlib.h>

// This line is only for the function CdinfoReceive. 
#include "MeshStruc.h" 
// Belows are for SparseView
#include "src/mat/impls/aij/seq/aij.h"

/*
VecReceive -- Receive a vector from Matlab 

Input parameters:
   PetscViewer -- a socket viewer obtained by PetscViewerSocketOpen
   VecType  -- can be VECMPI,VECSEQ or PETSC_NULL
Output parameters: 
   Vec -- The vector

Notes: 
   To send a vector from Matlab to Petsc, use the function send(portnum,A),
where portnum is obtained by the function openport().    

*/


// To use this function: 
// extern const int VecReceive(PetscViewer, VecType,Vec*);

int VecReceive(PetscViewer viewer, VecType outtype,Vec *newvec)
{
 
     MPI_Comm       comm;
     PetscMPIInt    rank,size,l=10,hostrank = 0,localsize;
     PetscInt       Istart,Iend,i;
     int            fd;    
     PetscErrorCode ierr; 
     Vec            vec;
     PetscScalar    *r,*localvec;
     PetscTruth     sameMPI,sameNULL;
     MPI_Status     status;
     
     // Get rank from viewer
     PetscObjectGetComm((PetscObject)viewer,&comm);
     MPI_Comm_rank(comm,&rank);
     MPI_Comm_size(comm,&size);

   
    if (rank==hostrank) {
     PetscInt        m,n;
 
     PetscViewerBinaryGetDescriptor(viewer,&fd); 
     // The matrix size is m*n
     ierr = PetscBinaryRead(fd,&m,1,PETSC_INT);CHKERRQ(ierr);
     ierr = PetscBinaryRead(fd,&n,1,PETSC_INT);CHKERRQ(ierr);
     l    = m*n;
     r    = malloc(l*sizeof(PetscScalar));
 
     // Followed by the data
     ierr = PetscBinaryRead(fd,r,l,PETSC_SCALAR);CHKERRQ(ierr);

     }

     // Broadcast the size of the vec    
     MPI_Bcast(&l,1,MPI_INT,hostrank,comm);

     ierr = PetscStrcmp(outtype,VECMPI,&sameMPI);CHKERRQ(ierr);
     ierr = PetscStrcmp(outtype,PETSC_NULL,&sameNULL);CHKERRQ(ierr);

     if (sameMPI||((l>1)&&(sameNULL))){

     	// Now create the MPI vec
        ierr = VecCreate(PETSC_COMM_WORLD,&vec);CHKERRQ(ierr);
     	ierr = VecSetSizes(vec,PETSC_DECIDE,l);CHKERRQ(ierr);
     	ierr = VecSetType(vec, VECMPI);CHKERRQ(ierr); 
        ierr = VecGetOwnershipRange(vec,&Istart,&Iend);CHKERRQ(ierr); 
        localsize = Iend-Istart;
        localvec  = malloc(localsize*sizeof(PetscScalar));
      
        if (rank!=hostrank){ 
             MPI_Send(&Istart,1,MPI_INT,hostrank,10*rank+0,comm);
             MPI_Send(&Iend  ,1,MPI_INT,hostrank,10*rank+1,comm);
             MPI_Recv(localvec,localsize,MPI_DOUBLE,hostrank,10*rank+2,comm,&status);

         }

        if (rank==hostrank){
              PetscInt localIstart,localIend;
              for (i=0;i<size;i++){
                 if (i!=hostrank){
                    MPI_Recv(&localIstart,1,MPI_INT,i,10*i+0,comm,&status);
                    MPI_Recv(&localIend  ,1,MPI_INT,i,10*i+1,comm,&status);
                    MPI_Send(r+localIstart,localIend-localIstart ,MPI_DOUBLE,i,10*i+2,comm);
                 }else{
                     PetscMemcpy(localvec,r,localsize*sizeof(PetscScalar));             
                     }
              }
             free(r);
         }

        PetscInt ix[localsize],i;
       	for(i=0;i<localsize;i++){ix[i]=Istart+i; }
        VecSetValues(vec,localsize,ix,localvec,INSERT_VALUES);
     	VecAssemblyBegin(vec);
     	VecAssemblyEnd(vec);
     
        
        //free(ix);
        free(localvec);
        *newvec = vec; 
         
     }else{
     	// Now we want to create a SEQ vector in 
     	// every processor, and all have the same values.         
        if(rank!=hostrank){ r    = malloc(l*sizeof(PetscScalar));}
     	MPI_Bcast(r,l,MPI_DOUBLE,hostrank,comm);
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,l,r,newvec);CHKERRQ(ierr);
             
           }
return(0);
}



///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*
MatReceiveTranspose -- Receive a sparse matrix from Matlab 

Input parameters:
   PetscViewer -- a socket viewer obtained by PetscViewerSocketOpen
   MatType  -- Has to be MATMPIAIJ
Output parameters: 
   Mat -- The matrix

Notes: 
   To send a matrix from Matlab to Petsc, use the function send(portnum,A),
where portnum is obtained by the function openport() and A HAS TO BE A SPARSE MATRIX.
Also, because Matlab stores sparse matrices column-oriented while Petsc does it 
row-oriented (in general). The matrix that Petscs receives is actually transpose(A). 
This is more efficient for matrix-vector multiplications.    

*/
// To use this function: 
//extern const int MatReceiveTranspose(PetscViewer, MatType, Mat*);

int MatReceiveTranspose(PetscViewer viewer, MatType outtype,Mat *newmat)
{

     MPI_Comm       comm;
     PetscMPIInt    rank,size,hostrank = 0,m,n;
     PetscInt       *Jc,*Ir,i,nnz,Istart,Iend,*localJc,*localIr,localnnz,Jc0;;
     int            fd;    
     PetscErrorCode ierr; 
     Mat            mat,matlocal;
     PetscScalar    *Pr,*localPr;
     MPI_Status     status;
      
     // Get rank from viewer
     PetscObjectGetComm((PetscObject)viewer,&comm);
     MPI_Comm_rank(comm,&rank);
     MPI_Comm_size(comm,&size);
   
  
    if (rank==hostrank) {
     
 
     PetscViewerBinaryGetDescriptor(viewer,&fd); 
     // The matrix size is m*n
     ierr = PetscBinaryRead(fd,&m,1,PETSC_INT);CHKERRQ(ierr);
     ierr = PetscBinaryRead(fd,&n,1,PETSC_INT);CHKERRQ(ierr);
     ierr = PetscBinaryRead(fd,&nnz,1,PETSC_INT);CHKERRQ(ierr);
    
     Jc    = malloc((n+1)*sizeof(PetscInt));
     Ir    = malloc(nnz*sizeof(PetscInt));
     Pr    = malloc(nnz*sizeof(PetscScalar));

     ierr = PetscBinaryRead(fd,Jc,n+1,PETSC_INT);CHKERRQ(ierr);
     ierr = PetscBinaryRead(fd,Ir,nnz,PETSC_INT);CHKERRQ(ierr);
     ierr = PetscBinaryRead(fd,Pr,nnz,PETSC_SCALAR);CHKERRQ(ierr);

     }
    
     // Broadcast the size of the mat
     MPI_Bcast(&m,1,MPI_INT,hostrank,comm);
     MPI_Bcast(&n,1,MPI_INT,hostrank,comm);
     MPI_Bcast(&nnz,1,MPI_INT,hostrank,comm);
     // Create a MPI matrix and get the row ownership range
     ierr = MatCreate(PETSC_COMM_WORLD,&mat);CHKERRQ(ierr);CHKERRQ(ierr);  
     ierr = MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,n,m);CHKERRQ(ierr);     ierr = MatSetType(mat,MATMPIAIJ);CHKERRQ(ierr);CHKERRQ(ierr);
     ierr = MatGetOwnershipRange(mat,&Istart,&Iend);CHKERRQ(ierr); 



   if (rank == hostrank){
       for (i=0;i<size;i++){
            PetscInt    localIstart,localIend,localnnz;          
            if(i!=hostrank){                  
            	MPI_Recv(&localIstart,1, MPI_INT, i, 10*i  ,comm,&status);
            	MPI_Recv(&localIend  ,1, MPI_INT, i, 10*i+1,comm,&status);
                localnnz = Jc[localIend]-Jc[localIstart] ;          
                MPI_Send(Jc+localIstart,localIend-localIstart+1,MPI_INT,i,10*i+2,comm);
                MPI_Send(Ir+Jc[localIstart],localnnz,MPI_INT,i,10*i+3,comm);
                MPI_Send(Pr+Jc[localIstart],localnnz,MPI_DOUBLE,i,10*i+4,comm);
             }       
         }       
     }

 
      localJc  = malloc((Iend-Istart+1)*sizeof(PetscInt));
      if (rank!=hostrank){
      	   MPI_Send(&Istart,1,MPI_INT,hostrank,10*rank,comm);
      	   MPI_Send(&Iend  ,1,MPI_INT,hostrank,10*rank+1,comm);
           MPI_Recv(localJc,Iend-Istart+1, MPI_INT,hostrank, 10*rank+2,comm,&status);

           localnnz = localJc[Iend-Istart]-localJc[0];
           localIr  = malloc(localnnz*sizeof(PetscInt));
           localPr  = malloc(localnnz*sizeof(PetscScalar));
      
      
           Jc0 = localJc[0];
           for (i=0;i<Iend-Istart+1;i++){localJc[i]=localJc[i]-Jc0;}
      	   MPI_Recv(localIr,localnnz, MPI_INT,hostrank, 10*rank+3,comm,&status);
      	   MPI_Recv(localPr,localnnz, MPI_DOUBLE,hostrank, 10*rank+4,comm,&status);  
        }else{

           PetscMemcpy(localJc,Jc,(Iend-Istart+1)*sizeof(PetscInt));
           localnnz = Jc[Iend-Istart]-Jc[0]; 
           localIr  = malloc(localnnz*sizeof(PetscInt));
           localPr  = malloc(localnnz*sizeof(PetscScalar));
           PetscMemcpy(localIr,Ir,localnnz*sizeof(PetscInt));
           PetscMemcpy(localPr,Pr,localnnz*sizeof(PetscScalar)); 
      
           free(Pr);
           free(Ir);
           free(Jc);

          }


   
      ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,Iend-Istart,m,localJc,localIr,localPr,&matlocal);CHKERRQ(ierr); 
      ierr = MatAssemblyBegin(matlocal,MAT_FINAL_ASSEMBLY);
      ierr = MatAssemblyEnd(matlocal,MAT_FINAL_ASSEMBLY);

      MatMerge(PETSC_COMM_WORLD,matlocal,PETSC_DECIDE,MAT_INITIAL_MATRIX,&mat);    
     *newmat = mat;
 
           free(localPr);
           free(localIr);
           free(localJc);
      


return(0);
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*
ScalarReceive -- Receive a scalar array from Matlab 

Input parameters:
   PetscViewer -- a socket viewer obtained by PetscViewerSocketOpen
   
Output parameters: 
   a -- The pointer to a scalar array

Notes: 
   To send a scalar from Matlab to Petsc, use the function send(portnum,a),
where portnum is obtained by the function openport().

*/


// To use this function: 
// extern const int ScalarReceive(PetscViewer, PetscScalar**);

int ScalarReceive(PetscViewer viewer, PetscScalar **a)
{
 
     MPI_Comm       comm;
     PetscMPIInt    rank,size,l=0,hostrank = 0;
     int            fd;    
     PetscErrorCode ierr; 
     PetscScalar    *r;
   
     
     // Get rank from viewer
     PetscObjectGetComm((PetscObject)viewer,&comm);
     MPI_Comm_rank(comm,&rank);
     MPI_Comm_size(comm,&size);

    
     if (rank==hostrank) {
         PetscInt        m,n;
 
         PetscViewerBinaryGetDescriptor(viewer,&fd); 
         // The matrix size is m*n
     	 ierr = PetscBinaryRead(fd,&m,1,PETSC_INT);CHKERRQ(ierr);
     	 ierr = PetscBinaryRead(fd,&n,1,PETSC_INT);CHKERRQ(ierr);
     	 l    = m*n;
     }
     
     MPI_Bcast(&l,1,MPI_INT,hostrank,comm);
     r    = malloc(l*sizeof(PetscScalar));

     if (rank==hostrank) {
         // Followed by the data
         ierr = PetscBinaryRead(fd,r,l,PETSC_SCALAR);CHKERRQ(ierr);
     }
     // Broadcast the data of the array
     MPI_Bcast(r,l,MPI_DOUBLE,hostrank,comm);
 
    
     *a  = r;
     

return(0);
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*
IntReceive -- Receive a integer array from Matlab 

Input parameters:
   PetscViewer -- a socket viewer obtained by PetscViewerSocketOpen
   
Output parameters: 
   a -- The pointer to an integer array

Notes: 
   To send an integer from Matlab to Petsc, use the function send(portnum,a),
where portnum is obtained by the function openport().

*/


// To use this function: 
// extern const int IntReceive(PetscViewer, PetscInt**);

int IntReceive(PetscViewer viewer, PetscInt **a)
{ 
 
     MPI_Comm       comm;
     PetscMPIInt    rank,size,l=0,hostrank = 0,i;
     int            fd;    
     PetscErrorCode ierr; 
     PetscScalar    *r;
     PetscInt       *ri;
   
     
     // Get rank from viewer
     PetscObjectGetComm((PetscObject)viewer,&comm);
     MPI_Comm_rank(comm,&rank);
     MPI_Comm_size(comm,&size);

    
     if (rank==hostrank) {
         PetscInt        m,n;
 
         PetscViewerBinaryGetDescriptor(viewer,&fd); 
         // The matrix size is m*n
         ierr = PetscBinaryRead(fd,&m,1,PETSC_INT);CHKERRQ(ierr);
       	 ierr = PetscBinaryRead(fd,&n,1,PETSC_INT);CHKERRQ(ierr);
         l    = m*n;
      }
     
      MPI_Bcast(&l,1,MPI_INT,hostrank,comm);
      r    = malloc(l*sizeof(PetscScalar));
      ri   = malloc(l*sizeof(PetscInt));
   
    
     if (rank==hostrank) {
        // Followed by the data
        ierr = PetscBinaryRead(fd,r,l,PETSC_SCALAR);CHKERRQ(ierr);
        // cast to integer. This is not very good!
        for (i=0;i<l;i++){*(ri+i)=(PetscInt)(*(r+i)); }
     }
  
     // Broadcast the data of the array
     MPI_Bcast(ri,l,MPI_INT,hostrank,comm); 

      
     *a  = ri;
     free(r);
     

return(0);
}

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//  Below are functions for my own use
// 
//

int CdinfoReceive(PetscViewer viewer, CoordInfo **cdinfoind)
{
  PetscErrorCode ierr;
  CoordInfo     *cdinfo;
  
  PetscInt     *dim;
  PetscInt  *l;
  PetscInt  *r;
  PetscInt     *n;
  PetscScalar  *gsize;
  PetscInt     *period;
 
  
  PetscInt    *gridspam;
  PetscInt    *nofgrid;

  cdinfo = malloc(sizeof(CoordInfo));
 
  ierr = IntReceive(viewer, &dim);CHKERRQ(ierr); 
  ierr = ScalarReceive(viewer, &l)  ;CHKERRQ(ierr); 
  ierr = ScalarReceive(viewer, &r)  ;CHKERRQ(ierr); 
  ierr = IntReceive(viewer, &n)  ;CHKERRQ(ierr); 
  ierr = ScalarReceive(viewer, &gsize);CHKERRQ(ierr); 
  ierr = IntReceive(viewer, &period);CHKERRQ(ierr); 
  ierr = IntReceive(viewer, &gridspam)  ;CHKERRQ(ierr); 
  ierr = IntReceive(viewer, &nofgrid);CHKERRQ(ierr); 


    cdinfo->dim      = dim;
    cdinfo->l        = (PetscInt3*)l;
    cdinfo->r        = (PetscInt3*)r;
    cdinfo->grid     = (PetscInt3*)n;
    cdinfo->gsize    = (PetscScalar3*)gsize;
    cdinfo->period   = (PetscInt3*)period;
    cdinfo->gridspam = (PetscInt33*)gridspam;
    cdinfo->nofgrid  = (PetscInt4*)nofgrid;
   

 
 
  *cdinfoind = cdinfo;
 return 0;
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*
  This function sends a Sparse MPI matrix to a socket. I write this because MatView
cannot send sparse mpi matrices to a socket. 
  In Matlab side, use SparseReceive to receive the sparse matrix.

Tzu-Chen Liang 2-14-2006

*/
int SparseView(Mat mat, PetscViewer viewer)
{
    PetscErrorCode ierr;
    PetscInt       Istart,Iend,Lsize;
    Mat            lmat;
    PetscInt       *ir,*jc,nnz;
    PetscScalar    *pr;
    Mat_SeqAIJ     *aij;

    ierr   = MatGetOwnershipRange(mat,&Istart,&Iend);CHKERRQ(ierr);  
    Lsize  = Iend - Istart;

 

    MatGetLocalMat(mat,MAT_INITIAL_MATRIX,&lmat); 
    aij = (Mat_SeqAIJ*)lmat->data;
    ir  = aij->i;
    jc  = aij->j;
    pr  = aij->a;
    nnz  =  *(ir+Lsize);

 

    PetscIntView(Lsize+1,ir,viewer); 
    PetscIntView(nnz,jc,viewer); 
    PetscScalarView(nnz,pr,viewer);

    ierr = MatDestroy(lmat);CHKERRQ(ierr);
  
return 0;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*
  This function sends a Sparse MPI matrix to 3 files.
  In Matlab side, use SparseLoad to receive the sparse matrix.
  

Tzu-Chen Liang 9-7-2007

*/
int SparseSave(Mat mat, char fname[])
{
    PetscErrorCode ierr;
    PetscInt       Istart,Iend,Lsize,rank,size,i,k;
    Mat            lmat;
    PetscInt       *ir,*jc,nnz;
    PetscScalar    *pr;
    Mat_SeqAIJ     *aij;
    FILE           *irfd,*jcfd,*prfd; 
    MPI_Status     status;
    PetscInt       *irtemp,*jctemp;
    PetscScalar    *prtemp;
    char           fnameir[50],fnamejc[50],fnamepr[50];

    ierr   = MatGetOwnershipRange(mat,&Istart,&Iend);CHKERRQ(ierr);  
    Lsize  = Iend - Istart;

    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


    MatGetLocalMat(mat,MAT_INITIAL_MATRIX,&lmat); 
    aij = (Mat_SeqAIJ*)lmat->data;
    ir  = aij->i;
    jc  = aij->j;
    pr  = aij->a;
    nnz  =  *(ir+Lsize);


    sprintf(fnameir,"%sirfiles",fname);
    sprintf(fnamejc,"%sjcfiles",fname);
    sprintf(fnamepr,"%sprfiles",fname);

  // Save the ir,jc,pr files  
    if(rank==0){

        PetscFOpen(PETSC_COMM_WORLD,fnameir,"w",&irfd);
        PetscFOpen(PETSC_COMM_WORLD,fnamejc,"w",&jcfd);
        PetscFOpen(PETSC_COMM_WORLD,fnamepr,"w",&prfd);
        //  PetscFOpen(PETSC_COMM_WORLD,"irfiles","w",&irfd);
        //  PetscFOpen(PETSC_COMM_WORLD,"jcfiles","w",&jcfd);
        //  PetscFOpen(PETSC_COMM_WORLD,"prfiles","w",&prfd);
    
     	  for(i=0;i<Lsize+1;i++){PetscFPrintf(PETSC_COMM_WORLD,irfd,"%d\n",*(ir+i));}
          for(i=0;i<nnz;i++){PetscFPrintf(PETSC_COMM_WORLD,jcfd,"%d\n",*(jc+i));}
          for(i=0;i<nnz;i++){PetscFPrintf(PETSC_COMM_WORLD,prfd,"%12.8f\n",*(pr+i));}
        
          for(k=1;k<size;k++){
  
             MPI_Recv(&Lsize,1,MPI_INT,k,10*k,PETSC_COMM_WORLD,&status);
             MPI_Recv(&nnz,1,MPI_INT,k,10*k+1,PETSC_COMM_WORLD,&status);  

             irtemp    = malloc((Lsize+1)*sizeof(PetscInt));
             jctemp    = malloc(nnz*sizeof(PetscInt));
             prtemp    = malloc(nnz*sizeof(PetscScalar));

      	     MPI_Recv(irtemp,Lsize+1,MPI_INT   ,k,10*k+2,PETSC_COMM_WORLD,&status);
      	     MPI_Recv(jctemp,nnz    ,MPI_INT   ,k,10*k+3,PETSC_COMM_WORLD,&status);
      	     MPI_Recv(prtemp,nnz    ,MPI_DOUBLE,k,10*k+4,PETSC_COMM_WORLD,&status);

        
     	     for(i=0;i<Lsize+1;i++){PetscFPrintf(PETSC_COMM_WORLD,irfd,"%d\n",*(irtemp+i));}
             for(i=0;i<nnz;i++){PetscFPrintf(PETSC_COMM_WORLD,jcfd,"%d\n",*(jctemp+i));}
             for(i=0;i<nnz;i++){PetscFPrintf(PETSC_COMM_WORLD,prfd,"%12.8f\n",*(prtemp+i));}
 
             free(irtemp);
             free(jctemp);
             free(prtemp);
   

           }
          PetscFClose(PETSC_COMM_WORLD,irfd);
          PetscFClose(PETSC_COMM_WORLD,jcfd);
          PetscFClose(PETSC_COMM_WORLD,prfd);
 
   }else{
      
             MPI_Send(&Lsize,1      ,MPI_INT   ,0,rank*10  ,PETSC_COMM_WORLD);
             MPI_Send(&nnz  ,1      ,MPI_INT   ,0,rank*10+1,PETSC_COMM_WORLD);
   	     MPI_Send(ir    ,Lsize+1,MPI_INT   ,0,rank*10+2,PETSC_COMM_WORLD);
    	     MPI_Send(jc    ,nnz    ,MPI_INT   ,0,rank*10+3,PETSC_COMM_WORLD);
    	     MPI_Send(pr    ,nnz    ,MPI_DOUBLE,0,rank*10+4,PETSC_COMM_WORLD);

   }

   // PetscIntView(Lsize+1,ir,viewer); 
   // PetscIntView(nnz,jc,viewer); 
   // PetscScalarView(nnz,pr,viewer);
   

    ierr = MatDestroy(lmat);CHKERRQ(ierr);
  
return 0;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*
  This function sends a vec to a file (Ascii)
  In Matlab side, use VecLoad to load the vector.


Tzu-Chen Liang 9-10-2007

*/
int VecSave(Vec vec, char fname[])
{
    PetscErrorCode ierr;
    PetscInt       Istart,Iend,Lsize,rank,size,i,k;
  
    FILE           *vecfd; 
    MPI_Status     status;
    PetscScalar    *vectemp, *vechead;

    ierr   = VecGetOwnershipRange(vec,&Istart,&Iend);CHKERRQ(ierr);  
    Lsize  = Iend - Istart;

    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);


    ierr   = VecGetArray(vec,&vechead);CHKERRQ(ierr); 



  // Save the file  
    if(rank==0){


          PetscFOpen(PETSC_COMM_WORLD,fname,"w",&vecfd);
          
     
          for(i=0;i<Lsize;i++){PetscFPrintf(PETSC_COMM_WORLD,vecfd,"%12.8f\n",*(vechead+i));}
        
          for(k=1;k<size;k++){
  
             MPI_Recv(&Lsize,1,MPI_INT,k,10*k,PETSC_COMM_WORLD,&status);
             vectemp    = malloc(Lsize*sizeof(PetscScalar));
      	     MPI_Recv(vectemp,Lsize    ,MPI_DOUBLE,k,10*k+1,PETSC_COMM_WORLD,&status);
    
             for(i=0;i<Lsize;i++){PetscFPrintf(PETSC_COMM_WORLD,vecfd,"%12.8f\n",*(vectemp+i));}
 
             free(vectemp);


           }
          PetscFClose(PETSC_COMM_WORLD,vecfd);
 
 
   }else{
      
             MPI_Send(&Lsize  ,1         ,MPI_INT   ,0,rank*10  ,PETSC_COMM_WORLD);
             MPI_Send(vechead ,Lsize     ,MPI_DOUBLE,0,rank*10+1,PETSC_COMM_WORLD);

   }

  
  
return 0;
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*
  This function sends a vec to a binary file 
  and can be loaded later by the function Vecloadfromfile
  the data format is : m,n,vecdata...

Tzu-Chen Liang 9-11-2007

*/
int VecSavebin(Vec vec, char fname[])
{
    PetscErrorCode ierr;
    PetscInt       Istart,Iend,Lsize,rank,size,i,k;
    PetscInt       m,n;
  
    FILE           *vecfd; 
    MPI_Status     status;
    PetscScalar    *vectemp, *vechead;

    ierr   = VecGetOwnershipRange(vec,&Istart,&Iend);CHKERRQ(ierr);  
    Lsize  = Iend - Istart;

    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

    VecGetSize(vec,&m);
    n  = 1;

    ierr   = VecGetArray(vec,&vechead);CHKERRQ(ierr); 



  // Save the file  
    if(rank==0){

          PetscBinaryOpen(fname,FILE_MODE_WRITE,&vecfd);
   
         //first save the size
         PetscBinaryWrite(vecfd,&m,1,PETSC_INT,PETSC_FALSE);        
         PetscBinaryWrite(vecfd,&n,1,PETSC_INT,PETSC_FALSE);        


         //followed by the data  
         PetscBinaryWrite(vecfd,vechead,Lsize,PETSC_SCALAR,PETSC_FALSE);
          
        
          for(k=1;k<size;k++){
  
             MPI_Recv(&Lsize,1,MPI_INT,k,10*k,PETSC_COMM_WORLD,&status);
             vectemp    = malloc(Lsize*sizeof(PetscScalar));
      	     MPI_Recv(vectemp,Lsize    ,MPI_DOUBLE,k,10*k+1,PETSC_COMM_WORLD,&status);
    
             PetscBinaryWrite(vecfd,vechead,Lsize,PETSC_SCALAR,PETSC_FALSE);
             
 
             free(vectemp);


           }
          
          PetscBinaryClose(vecfd);
 
   }else{
      
             MPI_Send(&Lsize  ,1         ,MPI_INT   ,0,rank*10  ,PETSC_COMM_WORLD);
             MPI_Send(vechead ,Lsize     ,MPI_DOUBLE,0,rank*10+1,PETSC_COMM_WORLD);

   }

  
  
return 0;
}



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/* 
  This function get the local ir,jc and pr data from an MPI matrix. 
  We can access the data but should not change it!

  Tzu-Chen Liang 2-19-2006
*/
int MatGetDataLocal(Mat mat,PetscInt **Ir, PetscInt **Jc,PetscScalar **Pr, PetscInt* nnz)
{
  PetscErrorCode ierr;
  Mat            lmat;
  PetscInt       Istart,Iend,Lsize;
  Mat_SeqAIJ     *aij;

  ierr   = MatGetOwnershipRange(mat,&Istart,&Iend);CHKERRQ(ierr);  
  Lsize  = Iend - Istart;
 

  MatGetLocalMat(mat,MAT_INITIAL_MATRIX,&lmat); 
  aij = (Mat_SeqAIJ*)lmat->data;
  *Ir  = aij->i;
  *Jc  = aij->j;
  *Pr  = aij->a;
  nnz  =  (*Ir+Lsize);


 return 0;
}
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/*
VecLoadfromFile -- Load a vector which is written by Matlab 

Input parameters:
   finemame -- the file name!
   VecType  -- can be VECMPI,VECSEQ or PETSC_NULL
Output parameters: 
   Vec -- The vector

Notes: 
   To send a vector from Matlab to Petsc, use the function send2file(filename,A),
    

*/


// To use this function: 
// extern const int VecReceive(char[], VecType,Vec*);

int VecLoadfromFile(MPI_Comm comm,char filename[], VecType outtype,Vec *newvec)
{
 
    
     PetscMPIInt    rank,size,l=10,hostrank = 0,localsize;
     PetscInt       Istart,Iend,i;
     int            fd;    
     PetscErrorCode ierr; 
     Vec            vec;
     PetscScalar    *r,*localvec;
     PetscTruth     sameMPI,sameNULL;
     MPI_Status     status;
     
     // Get rank from viewer
     //PetscObjectGetComm((PetscObject)viewer,&comm);
     MPI_Comm_rank(comm,&rank);
     MPI_Comm_size(comm,&size);

   
    if (rank==hostrank) {
     PetscInt        m,n;
     PetscBinaryOpen(filename,FILE_MODE_READ,&fd);
   
     // The matrix size is m*n
     ierr = PetscBinaryRead(fd,&m,1,PETSC_INT);CHKERRQ(ierr);
     ierr = PetscBinaryRead(fd,&n,1,PETSC_INT);CHKERRQ(ierr);
     l    = m*n;
     r    = malloc(l*sizeof(PetscScalar));
 
     // Followed by the data
     ierr = PetscBinaryRead(fd,r,l,PETSC_SCALAR);CHKERRQ(ierr);

     ierr = PetscBinaryClose(fd);
     }

     // Broadcast the size of the vec    
     MPI_Bcast(&l,1,MPI_INT,hostrank,comm);

     ierr = PetscStrcmp(outtype,VECMPI,&sameMPI);CHKERRQ(ierr);
     ierr = PetscStrcmp(outtype,PETSC_NULL,&sameNULL);CHKERRQ(ierr);

     if (sameMPI||((l>1)&&(sameNULL))){

     	// Now create the MPI vec
        ierr = VecCreate(PETSC_COMM_WORLD,&vec);CHKERRQ(ierr);
     	ierr = VecSetSizes(vec,PETSC_DECIDE,l);CHKERRQ(ierr);
     	ierr = VecSetType(vec, VECMPI);CHKERRQ(ierr); 
        ierr = VecGetOwnershipRange(vec,&Istart,&Iend);CHKERRQ(ierr); 
        localsize = Iend-Istart;
        localvec  = malloc(localsize*sizeof(PetscScalar));
      
        if (rank!=hostrank){ 
             MPI_Send(&Istart,1,MPI_INT,hostrank,10*rank+0,comm);
             MPI_Send(&Iend  ,1,MPI_INT,hostrank,10*rank+1,comm);
             MPI_Recv(localvec,localsize,MPI_DOUBLE,hostrank,10*rank+2,comm,&status);

         }

        if (rank==hostrank){
              PetscInt localIstart,localIend;
              for (i=0;i<size;i++){
                 if (i!=hostrank){
                    MPI_Recv(&localIstart,1,MPI_INT,i,10*i+0,comm,&status);
                    MPI_Recv(&localIend  ,1,MPI_INT,i,10*i+1,comm,&status);
                    MPI_Send(r+localIstart,localIend-localIstart ,MPI_DOUBLE,i,10*i+2,comm);
                 }else{
                     PetscMemcpy(localvec,r,localsize*sizeof(PetscScalar));             
                     }
              }
             free(r);
         }

        PetscInt ix[localsize],i;
       	for(i=0;i<localsize;i++){ix[i]=Istart+i; }
        VecSetValues(vec,localsize,ix,localvec,INSERT_VALUES);
     	VecAssemblyBegin(vec);
     	VecAssemblyEnd(vec);
     
        
        //free(ix);
        free(localvec);
        *newvec = vec; 
         
     }else{
     	// Now we want to create a SEQ vector in 
     	// every processor, and all have the same values.         
        if(rank!=hostrank){ r    = malloc(l*sizeof(PetscScalar));}
     	MPI_Bcast(r,l,MPI_DOUBLE,hostrank,comm);
        ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,l,r,newvec);CHKERRQ(ierr);
             
           }
return(0);
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*
MatLoadTransposefromFile -- Load a sparse matrix from generated and stored by Matlab 

Input parameters:
   filename -- the file name
   MatType  -- Has to be MATMPIAIJ
Output parameters: 
   Mat -- The matrix

Notes: 
   To send a matrix from Matlab to Petsc, use the function send2file(filename,A),
and A HAS TO BE A SPARSE MATRIX. Also, because Matlab stores sparse matrices 
column-oriented while Petsc does it row-oriented (in general). The matrix that 
Petscs receives is actually transpose(A). This is more efficient for matrix-vector 
multiplications.    

*/
// To use this function: 
//extern const int MatLoadTranspose(char[], MatType, Mat*);

int MatLoadTransposefromFile(MPI_Comm comm,char filename[], MatType outtype,Mat *newmat)
{

    
     PetscMPIInt    rank,size,hostrank = 0,m,n;
     PetscInt       *Jc,*Ir,i,nnz,Istart,Iend,*localJc,*localIr,localnnz,Jc0;;
     int            fd;    
     PetscErrorCode ierr; 
     Mat            mat,matlocal;
     PetscScalar    *Pr,*localPr;
     MPI_Status     status;
      
     // Get rank from viewer
     //PetscObjectGetComm((PetscObject)viewer,&comm);
     MPI_Comm_rank(comm,&rank);
     MPI_Comm_size(comm,&size);
   
  
    if (rank==hostrank) {
     
     ierr = PetscBinaryOpen(filename,FILE_MODE_READ,&fd);
      
     // The matrix size is m*n
     ierr = PetscBinaryRead(fd,&m,1,PETSC_INT);CHKERRQ(ierr);
     ierr = PetscBinaryRead(fd,&n,1,PETSC_INT);CHKERRQ(ierr);
     ierr = PetscBinaryRead(fd,&nnz,1,PETSC_INT);CHKERRQ(ierr);
    
     Jc    = malloc((n+1)*sizeof(PetscInt));
     Ir    = malloc(nnz*sizeof(PetscInt));
     Pr    = malloc(nnz*sizeof(PetscScalar));

     ierr = PetscBinaryRead(fd,Jc,n+1,PETSC_INT);CHKERRQ(ierr);
     ierr = PetscBinaryRead(fd,Ir,nnz,PETSC_INT);CHKERRQ(ierr);
     ierr = PetscBinaryRead(fd,Pr,nnz,PETSC_SCALAR);CHKERRQ(ierr);
   
     ierr = PetscBinaryClose(fd);

     }
    
     // Broadcast the size of the mat
     MPI_Bcast(&m,1,MPI_INT,hostrank,comm);
     MPI_Bcast(&n,1,MPI_INT,hostrank,comm);
     MPI_Bcast(&nnz,1,MPI_INT,hostrank,comm);
     // Create a MPI matrix and get the row ownership range
     ierr = MatCreate(PETSC_COMM_WORLD,&mat);CHKERRQ(ierr);CHKERRQ(ierr);  
     ierr = MatSetSizes(mat,PETSC_DECIDE,PETSC_DECIDE,n,m);CHKERRQ(ierr);     ierr = MatSetType(mat,MATMPIAIJ);CHKERRQ(ierr);CHKERRQ(ierr);
     ierr = MatGetOwnershipRange(mat,&Istart,&Iend);CHKERRQ(ierr); 



   if (rank == hostrank){
       for (i=0;i<size;i++){
            PetscInt    localIstart,localIend,localnnz;          
            if(i!=hostrank){                  
            	MPI_Recv(&localIstart,1, MPI_INT, i, 10*i  ,comm,&status);
            	MPI_Recv(&localIend  ,1, MPI_INT, i, 10*i+1,comm,&status);
                localnnz = Jc[localIend]-Jc[localIstart] ;          
                MPI_Send(Jc+localIstart,localIend-localIstart+1,MPI_INT,i,10*i+2,comm);
                MPI_Send(Ir+Jc[localIstart],localnnz,MPI_INT,i,10*i+3,comm);
                MPI_Send(Pr+Jc[localIstart],localnnz,MPI_DOUBLE,i,10*i+4,comm);
             }       
         }       
     }

 
      localJc  = malloc((Iend-Istart+1)*sizeof(PetscInt));
      if (rank!=hostrank){
      	   MPI_Send(&Istart,1,MPI_INT,hostrank,10*rank,comm);
      	   MPI_Send(&Iend  ,1,MPI_INT,hostrank,10*rank+1,comm);
           MPI_Recv(localJc,Iend-Istart+1, MPI_INT,hostrank, 10*rank+2,comm,&status);

           localnnz = localJc[Iend-Istart]-localJc[0];
           localIr  = malloc(localnnz*sizeof(PetscInt));
           localPr  = malloc(localnnz*sizeof(PetscScalar));
      
      
           Jc0 = localJc[0];
           for (i=0;i<Iend-Istart+1;i++){localJc[i]=localJc[i]-Jc0;}
      	   MPI_Recv(localIr,localnnz, MPI_INT,hostrank, 10*rank+3,comm,&status);
      	   MPI_Recv(localPr,localnnz, MPI_DOUBLE,hostrank, 10*rank+4,comm,&status);  
        }else{

           PetscMemcpy(localJc,Jc,(Iend-Istart+1)*sizeof(PetscInt));
           localnnz = Jc[Iend-Istart]-Jc[0]; 
           localIr  = malloc(localnnz*sizeof(PetscInt));
           localPr  = malloc(localnnz*sizeof(PetscScalar));
           PetscMemcpy(localIr,Ir,localnnz*sizeof(PetscInt));
           PetscMemcpy(localPr,Pr,localnnz*sizeof(PetscScalar)); 
      
           free(Pr);
           free(Ir);
           free(Jc);

          }


   
      ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,Iend-Istart,m,localJc,localIr,localPr,&matlocal);CHKERRQ(ierr); 
      ierr = MatAssemblyBegin(matlocal,MAT_FINAL_ASSEMBLY);
      ierr = MatAssemblyEnd(matlocal,MAT_FINAL_ASSEMBLY);

      MatMerge(PETSC_COMM_WORLD,matlocal,PETSC_DECIDE,MAT_INITIAL_MATRIX,&mat);    
     *newmat = mat;
 
           free(localPr);
           free(localIr);
           free(localJc);
      


return(0);
}


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
/*
IntLoadfromFile -- Load a integer array stored by Matlab 

Input parameters:
   filename -- the file name
   
Output parameters: 
   a -- The pointer to an integer array

Notes: 
   To send an integer from Matlab to Petsc, use the function send2file(filename,a),


*/


// To use this function: 
// extern const int IntLoadfromFile(PETSC_COMM_WORLD, filename, PetscInt**);
int IntLoadfromFile(MPI_Comm comm, char filename[], PetscInt **a)
{ 
 
     PetscMPIInt    rank,size,l=0,hostrank = 0,i;
     int            fd;    
     PetscErrorCode ierr; 
     PetscScalar    *r;
     PetscInt       *ri;
   

     MPI_Comm_rank(comm,&rank);
     MPI_Comm_size(comm,&size);

    
     if (rank==hostrank) {
         PetscInt        m,n;
         ierr = PetscBinaryOpen(filename,FILE_MODE_READ,&fd);
         
         // The matrix size is m*n
         ierr = PetscBinaryRead(fd,&m,1,PETSC_INT);CHKERRQ(ierr);
       	 ierr = PetscBinaryRead(fd,&n,1,PETSC_INT);CHKERRQ(ierr);
         l    = m*n;
         
      }
     
      MPI_Bcast(&l,1,MPI_INT,hostrank,comm);
      r    = malloc(l*sizeof(PetscScalar));
      ri   = malloc(l*sizeof(PetscInt));
   
    
     if (rank==hostrank) {
        // Followed by the data
        ierr = PetscBinaryRead(fd,r,l,PETSC_SCALAR);CHKERRQ(ierr);
        // cast to integer. This is not very good!
        for (i=0;i<l;i++){*(ri+i)=(PetscInt)(*(r+i)); }
        ierr = PetscBinaryClose(fd);
     }
  
     // Broadcast the data of the array
     MPI_Bcast(ri,l,MPI_INT,hostrank,comm); 

      
     *a  = ri;
     free(r);
     

return(0);
}
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*
ScalarLoadfromFile -- Load a scalar array stored by Matlab 

Input parameters:
   filename -- the filename (binary)
   
Output parameters: 
   a -- The pointer to a scalar array

Notes: 
   To save a scalar from Matlab to Petsc, use the function send2file(filename,a),


*/


// To use this function: 
// extern const int ScalarLoadfromFile(PETSC_COMM_WORLD,filename, PetscScalar**);

int ScalarLoadfromFile(MPI_Comm comm, char filename[], PetscScalar **a)
{
 
     PetscMPIInt    rank,size,l=0,hostrank = 0;
     int            fd;    
     PetscErrorCode ierr; 
     PetscScalar    *r;
   

     MPI_Comm_rank(comm,&rank);
     MPI_Comm_size(comm,&size);

    
     if (rank==hostrank) {
         PetscInt        m,n;
 
         ierr = PetscBinaryOpen(filename,FILE_MODE_READ,&fd);
         // The matrix size is m*n
     	 ierr = PetscBinaryRead(fd,&m,1,PETSC_INT);CHKERRQ(ierr);
     	 ierr = PetscBinaryRead(fd,&n,1,PETSC_INT);CHKERRQ(ierr);
     	 l    = m*n;
     }
     
     MPI_Bcast(&l,1,MPI_INT,hostrank,comm);
     r    = malloc(l*sizeof(PetscScalar));

     if (rank==hostrank) {
         // Followed by the data
         ierr = PetscBinaryRead(fd,r,l,PETSC_SCALAR);CHKERRQ(ierr);
         ierr = PetscBinaryClose(fd);
     }
     // Broadcast the data of the array
     MPI_Bcast(r,l,MPI_DOUBLE,hostrank,comm);
 
    
     *a  = r;
     

return(0);
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

int CdinfoLoad(MPI_Comm comm,char cdinfofile[], CoordInfo **cdinfoind)
{


  PetscErrorCode ierr;
  CoordInfo     *cdinfo;
  
  PetscInt     *dim;
  PetscInt     *l;
  PetscInt     *r;
  PetscInt     *n;
  PetscScalar  *gsize;
  PetscInt     *period;
  PetscInt     *temp,i;
  PetscScalar  *temp2;
  
  PetscInt    *gridspam;
  PetscInt    *nofgrid;

  cdinfo = malloc(sizeof(CoordInfo));
  temp  = malloc(29*sizeof(PetscInt));
  gsize = malloc(3*sizeof(PetscScalar));
  

  ierr = IntLoadfromFile(comm,cdinfofile, &temp);CHKERRQ(ierr); 
 
   dim      = temp;
   l        = temp+1;
   r        = temp+4;
   n        = temp+7;
   period   = temp+10;
   gridspam = temp+13;
   nofgrid  = temp+25;



   *(gsize+0) = (*r-*l)*1.0/(*n); 
   *(gsize+1) = (*(r+1)-*(l+1))*1.0/(*(n+1));
   *(gsize+2) = (*(r+2)-*(l+2))*1.0/(*(n+2));

    cdinfo->dim      = dim;
    cdinfo->l        = (PetscInt3*)l;
    cdinfo->r        = (PetscInt3*)r;
    cdinfo->grid     = (PetscInt3*)n;
    cdinfo->gsize    = (PetscScalar3*)gsize;
    cdinfo->period   = (PetscInt3*)period;
    cdinfo->gridspam = (PetscInt33*)gridspam;
    cdinfo->nofgrid  = (PetscInt4*)nofgrid;
   


 
  *cdinfoind = cdinfo;
 return 0;
}
