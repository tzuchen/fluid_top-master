/*
  Here I store the iteration functions

  
  Tzu-Chen Liang 11-6-2006
*/

#include "petscksp.h"
#include "petscda.h"
#include "petsctime.h"
#include <math.h>
#include <stdlib.h>
#include "petscbt.h"
#include "IterationFunction.h"
#include "LargeVecFunction.h"
////////////////////////////////////////////////////////////////////////////////////////////////
/*
Perform backward averaging iteration 
 x    : data
 y    : next iteration

 pmax : number of points to be cached
 Istart, Iend: column  

 cacheInt    : (2*npt+1)*pmax*sizeof(PetscInt)
 cacheScalar : (  npt+1)*pmax*sizeof(PetscScalar)

Tzu-Chen Liang  11-25-2006

*/
int BackwardAverage(Vec *x, Vec *y, PetscInt *cacheInt, PetscScalar *cacheScalar, PetscInt n, PetscInt npt, PetscInt pmax, PetscInt Istart, PetscInt Iend,PetscScalar c){

 PetscErrorCode ierr;
 PetscInt       i, j, k=0, pi, pj, n2,m, puse, pgrid;  

 PetscInt       localsizex,localsizey, rowcount=0;
 PetscInt       *idp, *NzindJ, *NzindI;
 PetscScalar    dx,dy,dx2,dy2,CX,CY;
 PetscScalar    *pty,*pty0; 
 Vec            y0;
 PetscInt       nvec = 2;
 Vec            xvec[nvec];
 IS             ISfrom[nvec],ISto[nvec];
 VecScatter     ctxt[nvec];



 n2     = (PetscInt)(n*0.5);
 dx     = 1.0/n;
 dy     = 1.0/n;
 dx2    = dx/2-dx/1e6;
 dy2    = dy/2-dy/1e6;

  LargeVecCreate(x,nvec,xvec);
  NzindI = cacheInt;        //pmax
  NzindJ = cacheInt+pmax;   //pmax
  idp    = cacheInt;        //pmax
  pty0   = cacheScalar;     //pmax

  
  localsizex    = Iend-Istart;
  localsizey    = (PetscInt)(pmax*1.0/(localsizex+1))-3;
  if(localsizey>n2){localsizey =n2;}
 

  ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,(localsizex+1)*(localsizey+1),pty0,&y0);
  

while(rowcount<n2){

     if (n2-rowcount<=localsizey){localsizey =n2-rowcount;}      
     puse = localsizex*localsizey;
     pgrid = (localsizex+1)*(localsizey+1);
     k= 0;
     for(i=Istart;i<Iend+1;i++){
         for(j=rowcount;j<rowcount+localsizey+1;j++){
             CX = (PetscScalar)(i*dx);
             CY = (PetscScalar)(j*dy); 

              InverseStandardMap(&CX,&CY,c); 
             //InverseModifiedArnoldsCatMap(&CX,&CY);
 
             pi = (PetscInt)floor(CX*(PetscScalar)n);
             pj = (PetscInt)floor(CY*(PetscScalar)n);   
             if(pj>=n2) {SkewSymmetricPoint(&pi, &pj, n);}
             *(NzindI+k) = pi; 
             *(NzindJ+k) = pj;
             k++;    

          }
      }

    ierr =  VecDestroy(y0);CHKERRQ(ierr);
    ierr =  VecCreateSeqWithArray(PETSC_COMM_SELF,pgrid,pty0,&y0);CHKERRQ(ierr);
   
    ISCreateGeneralWithIJ(MPI_COMM_SELF,*x,xvec,nvec,n2,pgrid, NzindI, NzindJ,ISfrom, ISto);
    LargeVecScatterCreate(xvec,ISfrom,y0,ISto ,ctxt,nvec); 
    ISArrayDestroy(ISfrom,nvec);
    ISArrayDestroy(ISto,nvec);
    LargeVecScatterBeginEnd(xvec,y0,INSERT_VALUES,SCATTER_FORWARD,ctxt,nvec);
    VecScatterArrayDestroy(ctxt,nvec);


    ierr =  VecGetArray(y0,&pty0);CHKERRQ(ierr);
    ierr =  VecGetArray(*y,&pty);CHKERRQ(ierr);
    m    = 0;
      for(i=0;i<localsizex;i++){
           for(j=0;j<localsizey;j++){             
              *(pty+i*n2+j+rowcount) = (*(pty0+i*(localsizey+1)+j)+
                                      *(pty0+i*(localsizey+1)+j+1)+
                                      *(pty0+(i+1)*(localsizey+1)+j)+
                                      *(pty0+(i+1)*(localsizey+1)+j+1))/4;         
               m++; 
           }
      }
      VecRestoreArray(y0,&pty0);
      VecRestoreArray(*y,&pty);

     rowcount = rowcount + localsizey;


}

 ierr =  VecDestroy(y0);CHKERRQ(ierr);
 VecArrayDestroy(xvec,nvec);
 

return 0;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/*
A new Strategy can handle large size problem (more than 4G variables)

*/
int BackwardAverageRL(Vec *x, Vec *y, PetscInt *cacheInt, PetscScalar *cacheScalar, PetscInt n, PetscInt npt, PetscInt pmax, PetscInt Istart, PetscInt Iend,PetscScalar c){

 PetscInt       rank,size;
 PetscErrorCode ierr;
 PetscInt       i, j, k=0, pi, pj, n2,n4 ,m, puse, pgrid,lx;  
 PetscInt       localsizex,localsizey, rowcount=0;
 PetscInt       k1,k2,pgrid1,pgrid2;
 PetscInt       *idy,*idp, *NzindJ;
 PetscScalar    dx,dy,dx2,dy2,CX,CY;

 PetscScalar    *pty, *pty0; 
 IS             isx1,isx2,isy1,isy2;
 VecScatter     ctx1,ctx2;
 Vec            y0;
 Vec            x1,x2;
 PetscScalar    *ptx1,*ptx2;
 PetscInt       size1,size2,col1,col2;
 

 MPI_Comm_size(PETSC_COMM_WORLD,&size);
 MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

 n2     = (PetscInt)(n*0.5);
 n4     = (PetscInt)(n*0.25);
 dx     = 1.0/n;
 dy     = 1.0/n;
 dx2    = dx/2-dx/1e6;
 dy2    = dy/2-dy/1e6;


  NzindJ = cacheInt;    //pmax
  idp    = cacheInt;    //pmax
  idy    = cacheInt   + pmax; 

  pty0   = cacheScalar   ; //pmax

  
  localsizex    = Iend-Istart;
  localsizey    = (PetscInt)(pmax*1.0/(localsizex+1))-2;
  if(localsizey>n2){localsizey =n2;}
 
  
  ierr =  VecGetArray(*x,&ptx1);CHKERRQ(ierr);
  ptx2 = ptx1;

  if(rank< size*0.5){lx =  localsizex*n2;}else{lx =0;}
  VecCreateMPIWithArray(PETSC_COMM_WORLD,lx,PETSC_DETERMINE,ptx1,&x1);
  if(rank< size*0.5){lx =  0;}else{lx =  localsizex*n2; }
  VecCreateMPIWithArray(PETSC_COMM_WORLD,lx,PETSC_DETERMINE,ptx2,&x2);

  VecGetSize(x1,&size1);
  VecGetSize(x2,&size2);
  col1 = (PetscInt)(size1*1.0/n2);
  col2 = (PetscInt)(size2*1.0/n2);

  ierr =  VecGetArray(*y,&pty);CHKERRQ(ierr);
  ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,(localsizex+1)*(localsizey+1),pty0,&y0);


while(rowcount<n2){

     if (n2-rowcount<=localsizey){localsizey =n2-rowcount;}      
     puse = localsizex*localsizey;
     pgrid = (localsizex+1)*(localsizey+1);
     k= 0;
     k1=0;
     k2=0; 
     for(i=Istart;i<Iend+1;i++){
        for(j=rowcount;j<rowcount+localsizey+1;j++){
           
             CX = (PetscScalar)(i*dx);
             CY = (PetscScalar)(j*dy); 
             InverseStandardMap(&CX,&CY,c); 
             pi = (PetscInt)floor(CX*n);
             pj = (PetscInt)floor(CY*n);   
             if(pj>=n2) {SkewSymmetricPoint(&pi, &pj, n);}

             if(pi<col1){
                  *(NzindJ+k1) =  (PetscInt)(n2*pi +  pj);
                  *(idy+k1)   =  k;
                   k1++;
              }else{
                  *(NzindJ+pgrid-k2-1) =  (PetscInt)(n2*(pi-col1)+pj);
                  *(idy+pgrid-k2-1)   =  k;
                   k2++;
              } 
             k++;    
          }
      }

      pgrid1 = k1;
      pgrid2 = k2;


    ierr =  ISCreateGeneralWithArray(PETSC_COMM_SELF,pgrid1,NzindJ,&isx1);CHKERRQ(ierr);
    ierr =  ISCreateGeneralWithArray(PETSC_COMM_SELF,pgrid2,NzindJ+pgrid1,&isx2);CHKERRQ(ierr);  
 
    ierr =  ISCreateGeneralWithArray(PETSC_COMM_SELF,pgrid1,idy,&isy1);CHKERRQ(ierr);
    ierr =  ISCreateGeneralWithArray(PETSC_COMM_SELF,pgrid2,idy+pgrid1,&isy2);CHKERRQ(ierr);

  
    ierr =  VecDestroy(y0);CHKERRQ(ierr);
    ierr =  VecCreateSeqWithArray(PETSC_COMM_SELF,pgrid,pty0,&y0);CHKERRQ(ierr);



    ierr =  VecScatterCreate(x1,isx1,y0,isy1,&ctx1);CHKERRQ(ierr);
    ierr =  VecScatterCreate(x2,isx2,y0,isy2,&ctx2);CHKERRQ(ierr);

    ierr =  VecScatterBegin(x1,y0,INSERT_VALUES,SCATTER_FORWARD,ctx1);CHKERRQ(ierr);
    ierr =  VecScatterEnd(x1,y0,INSERT_VALUES,SCATTER_FORWARD,ctx1);CHKERRQ(ierr);
    ierr =  VecScatterBegin(x2,y0,INSERT_VALUES,SCATTER_FORWARD,ctx2);CHKERRQ(ierr);
    ierr =  VecScatterEnd(x2,y0,INSERT_VALUES,SCATTER_FORWARD,ctx2);CHKERRQ(ierr);

    ierr =  VecScatterDestroy(ctx1);
    ierr =  VecScatterDestroy(ctx2);
    ierr =  VecGetArray(y0,&pty0);CHKERRQ(ierr);

      m = 0;
      for(i=0;i<localsizex;i++){
           for(j=0;j<localsizey;j++){
              *(pty+i*n2+j+rowcount) = (*(pty0+i*(localsizey+1)+j)+
                                        *(pty0+i*(localsizey+1)+j+1)+
                                        *(pty0+(i+1)*(localsizey+1)+j)+
                                        *(pty0+(i+1)*(localsizey+1)+j+1))/4;         
               m++; 
           }
      }

     VecRestoreArray(y0,&pty0);
     VecRestoreArray(*y,&pty);
     rowcount = rowcount + localsizey;
}



     VecDestroy(x1);
     VecDestroy(x2);


return 0;
}
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
/*
Perform backward center iteration 
 x    : data
 y    : next iteration

 pmax : number of points to be cached
 Istart, Iend: column  

 cacheInt    : (2*npt+1)*pmax*sizeof(PetscInt)
 cacheScalar : (  npt+1)*pmax*sizeof(PetscScalar)

Tzu-Chen Liang  11-25-2006

*/
int BackwardCenter(Vec *x, Vec *y, PetscInt *cacheInt, PetscScalar *cacheScalar, PetscInt n, PetscInt npt, PetscInt pmax, PetscInt Istart, PetscInt Iend,PetscScalar c){

 PetscErrorCode ierr;
 PetscInt       i, j, k=0, pi, pj, n2,m, puse, pgrid;  

 PetscInt       localsizex,localsizey, rowcount=0;
 PetscInt       *idp, *NzindJ, *NzindI;
 PetscScalar    dx,dy,dx2,dy2,CX,CY;
 PetscScalar    *pty,*pty0; 
 Vec            y0;
 PetscInt       nvec = 2;
 Vec            xvec[nvec];
 IS             ISfrom[nvec],ISto[nvec];
 VecScatter     ctxt[nvec];



 n2     = (PetscInt)(n*0.5);
 dx     = 1.0/n;
 dy     = 1.0/n;
 dx2    = dx/2-dx/1e6;
 dy2    = dy/2-dy/1e6;

  LargeVecCreate(x,nvec,xvec);
  NzindI = cacheInt;        //pmax
  NzindJ = cacheInt+pmax;   //pmax
  idp    = cacheInt;        //pmax
  pty0   = cacheScalar;     //pmax

  
  localsizex    = Iend-Istart;
  localsizey    = (PetscInt)(pmax*1.0/(localsizex+1))-3;
  if(localsizey>n2){localsizey =n2;}
 

  ierr = VecCreateSeqWithArray(PETSC_COMM_SELF,(localsizex)*(localsizey),pty0,&y0);
  

while(rowcount<n2){

     if (n2-rowcount<=localsizey){localsizey =n2-rowcount;}      
     puse = localsizex*localsizey;
     pgrid = localsizex*localsizey;
     k= 0;
     for(i=Istart;i<Iend;i++){
         for(j=rowcount;j<rowcount+localsizey;j++){
             CX = (PetscScalar)(i*dx+dx/2);
             CY = (PetscScalar)(j*dy+dy/2); 
             InverseStandardMap(&CX,&CY,c); 
             pi = (PetscInt)floor(CX*(PetscScalar)n);
             pj = (PetscInt)floor(CY*(PetscScalar)n);   
             if(pj>=n2) {SkewSymmetricPoint(&pi, &pj, n);}
             *(NzindI+k) = pi; 
             *(NzindJ+k) = pj;
             k++;    

          }
      }

    ierr =  VecDestroy(y0);CHKERRQ(ierr);
    ierr =  VecCreateSeqWithArray(PETSC_COMM_SELF,pgrid,pty0,&y0);CHKERRQ(ierr);
   
    ISCreateGeneralWithIJ(MPI_COMM_SELF,*x,xvec,nvec,n2,pgrid, NzindI, NzindJ,ISfrom, ISto);
    LargeVecScatterCreate(xvec,ISfrom,y0,ISto ,ctxt,nvec); 
    ISArrayDestroy(ISfrom,nvec);
    ISArrayDestroy(ISto,nvec);
    LargeVecScatterBeginEnd(xvec,y0,INSERT_VALUES,SCATTER_FORWARD,ctxt,nvec);
    VecScatterArrayDestroy(ctxt,nvec);


    ierr =  VecGetArray(y0,&pty0);CHKERRQ(ierr);
    ierr =  VecGetArray(*y,&pty);CHKERRQ(ierr);
    m    = 0;
      for(i=0;i<localsizex;i++){
           for(j=0;j<localsizey;j++){             
              *(pty+i*n2+j+rowcount) = *(pty0+i*localsizey+j);                                               
               m++; 
           }
      }
      VecRestoreArray(y0,&pty0);
      VecRestoreArray(*y,&pty);

     rowcount = rowcount + localsizey;


}

 ierr =  VecDestroy(y0);CHKERRQ(ierr);
 VecArrayDestroy(xvec,nvec);
 

return 0;
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

/*
int ForwardInterpolation(Vec *x, Vec *y, PetscInt *cacheInt, PetscScalar *cacheScalar, PetscInt n, PetscInt npt, PetscInt pmax, PetscInt Istart, PetscInt Iend){


 


   }

return 0;

}
*/
//////////////////////////////////////////////////////////////////////////////////////
int Smoothing(Vec *x, Vec *y, PetscScalar *cacheScalar, PetscInt *cacheInt , VecScatter *ctx,PetscInt n,DA myDA, PetscInt Istart, PetscInt Iend){


PetscScalar    C0,C1,C2;
PetscErrorCode ierr;
PetscInt       localsizex;
PetscInt       n2,i, j, k; 
Vec            lvecx,lvecy;
PetscScalar    **lvecptx,**lvecpty;
Vec            bcvec;
IS             isbc;
//VecScatter     ctx;
PetscScalar    *bcpt1,*bcpt2,*bcpt3,*bcpt4;



localsizex    = Iend-Istart;

  C0 = 1.0/8;
  C1 = 1*1.0/4;
  C2 = 1*3.0/16;

n2     = (PetscInt)(n*0.5);


DACreateLocalVector(myDA,&lvecx);
DACreateLocalVector(myDA,&lvecy);
DAGlobalToLocalBegin(myDA,*x,INSERT_VALUES,lvecx);
DAGlobalToLocalEnd(myDA,*x,INSERT_VALUES,lvecx);



VecGetArray2d(lvecx,localsizex+4,n2+4,0,0,&lvecptx);
VecGetArray2d(lvecy,localsizex+4,n2+4,0,0,&lvecpty);

// X direction smoothing
for(j=0;j<n2+4;j++){
   for(i=2;i<localsizex+2;i++){

   lvecpty[i][j] = C0*lvecptx[i][j]
                  +C1*lvecptx[i-1][j] +C1*lvecptx[i+1][j]
                  +C2*lvecptx[i-2][j]+C2*lvecptx[i+2][j];             
   }   
}


VecCreateSeqWithArray(MPI_COMM_SELF,localsizex*4,cacheScalar,&bcvec);
VecScatterBegin(*x,bcvec,INSERT_VALUES,SCATTER_FORWARD,*ctx);
VecScatterEnd(*x,bcvec,INSERT_VALUES,SCATTER_FORWARD,*ctx);


bcpt1 = cacheScalar;
bcpt2 = cacheScalar+localsizex;
bcpt3 = cacheScalar+localsizex*2;
bcpt4 = cacheScalar+localsizex*3;

k= 0;
for(i=2;i<localsizex+2;i++){
 
  lvecpty[i][0]= *(bcpt3+k);
  lvecpty[i][1]= *(bcpt4+k);
  lvecpty[i][n2+3]= *(bcpt2+k);
  lvecpty[i][n2+2]= *(bcpt1+k);

  k++;
}

// Y direction smoothing
for(j=2;j<n2+2;j++){
   for(i=2;i<localsizex+2;i++){
   lvecptx[i][j] = C0*lvecpty[i][j]
                  +C1*lvecpty[i][j-1] +C1*lvecpty[i][j+1]
                  +C2*lvecpty[i][j-2] +C2*lvecpty[i][j+2];               
   }   
}


VecRestoreArray2d(lvecy,localsizex+4,n2+4,0,0,&lvecpty);
VecRestoreArray2d(lvecx,localsizex+4,n2+4,0,0,&lvecptx);

DALocalToGlobal(myDA,lvecx,INSERT_VALUES,*x);

VecDestroy(bcvec);
VecDestroy(lvecx);
VecDestroy(lvecy);

return 0;
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
int SmoothingRL(Vec *x, Vec *y, PetscScalar *cacheScalar, PetscInt *cacheInt , VecScatter *ctx,PetscInt n, PetscInt Istart, PetscInt Iend){

PetscInt       rank,size;
PetscScalar    C0,C1,C2;
PetscErrorCode ierr;
PetscInt       n2,i, j, k, localsizex;
PetscInt       nvec; 

Vec            xbc,ybc,lvecx,lvecy;
PetscScalar    *xbcpt,*ybcpt,**lvecptx,**lvecpty;
PetscScalar    **xbcpt2,**ybcpt2;
PetscInt       *bcISarray;
PetscInt       left,right;
PetscInt       *bcindI,*bcindJ;

nvec  = 2;
IS             ISfrom[nvec],ISto[nvec];
VecScatter     ctxt[nvec];
Vec            xvec[nvec];



MPI_Comm_size(PETSC_COMM_WORLD,&size);
MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
localsizex = Iend-Istart;

  C0 = 1.0/8;
  C1 = 1.0/4;
  C2 = 3.0/16;

n2     = (PetscInt)(n*0.5);
xbcpt  = cacheScalar;
ybcpt  = cacheScalar;
bcindI = cacheInt;
bcindJ = cacheInt+ 4*n2;


LargeVecCreate(x,nvec,xvec);


if(rank==0){left = n-1;}
       else{left = Istart-1;}
if(rank==size-1){right = 0;}
       else{right = Iend;}


for(i=0;i<n2;i++){ *(bcindI+i)      = left-1;
                   *(bcindJ+i)      = i;        }
for(i=0;i<n2;i++){ *(bcindI+n2+i)   = left;
                   *(bcindJ+n2+i)   = i;        }
for(i=0;i<n2;i++){ *(bcindI+n2*2+i) = right;
                   *(bcindJ+n2*2+i) = i;        }
for(i=0;i<n2;i++){ *(bcindI+n2*3+i) = right+1;
                   *(bcindJ+n2*3+i) = i;        }


  VecCreateSeqWithArray(MPI_COMM_SELF,4*n2,cacheScalar,&xbc);
  ISCreateGeneralWithIJ(MPI_COMM_SELF,x,xvec,nvec,n2,4*n2, bcindI, bcindJ,ISfrom, ISto);
  LargeVecScatterCreate(xvec,ISfrom,xbc,ISto ,ctxt,nvec); 
  ISArrayDestroy(ISfrom,nvec);
  ISArrayDestroy(ISto,nvec);
  LargeVecScatterBeginEnd(xvec,xbc,INSERT_VALUES,SCATTER_FORWARD,ctxt,nvec);
  VecScatterArrayDestroy(ctxt,nvec);


  VecGetArray2d(*x,localsizex,n2,0,0,&lvecptx);
  VecGetArray2d(*y,localsizex,n2,0,0,&lvecpty);
  VecGetArray2d(xbc,4,n2,0,0,&xbcpt2);
  
  for(j=0;j<n2;j++){
     for(i=2;i<localsizex-2;i++){

     lvecpty[i][j] = C0*lvecptx[i][j]
                    +C1*lvecptx[i-1][j] +C1*lvecptx[i+1][j]
                    +C2*lvecptx[i-2][j]+C2*lvecptx[i+2][j];             
     }   
  }
 
  for(j=0;j<n2;j++){
  
      lvecpty[0][j] = C0*lvecptx[0][j]
                    +C1*xbcpt2[1][j] +C1*lvecptx[1][j]
                    +C2*xbcpt2[0][j]+C2*lvecptx[2][j];  
      lvecpty[1][j] = C0*lvecptx[1][j]
                    +C1*lvecptx[0][j] +C1*lvecptx[2][j]
                    +C2*xbcpt2[1][j]+C2*lvecptx[3][j];  

      lvecpty[localsizex-1][j] = C0*lvecptx[localsizex-1][j]
                                +C1*lvecptx[localsizex-2][j] +C1*xbcpt2[2][j]
                                +C2*lvecptx[localsizex-3][j] +C2*xbcpt2[3][j];   
      lvecpty[localsizex-2][j] = C0*lvecptx[localsizex-2][j]
                                +C1*lvecptx[localsizex-3][j] +C1*lvecptx[localsizex-1][j]
                                +C2*lvecptx[localsizex-4][j] +C2*xbcpt2[2][j];  

  }

  VecRestoreArray2d(xbc,4,n2,0,0,&xbcpt2);
  VecDestroy(xbc);
//////////////////////////////////////////////////////////////////////////////////
k = 0;
for(i=Istart;i<Iend;i++){ *(bcindI+k)      = n-i-1;
                          *(bcindJ+k)      = 1;           
                          k++;                      }
for(i=Istart;i<Iend;i++){ *(bcindI+k)      = n-i-1;
                          *(bcindJ+k)      = 0;
                          k++;                      }
for(i=Istart;i<Iend;i++){ *(bcindI+k)      = n-i-1;
                          *(bcindJ+k)      = n2-1;
                          k++;                      }
for(i=Istart;i<Iend;i++){ *(bcindI+k)      = n-i-1;
                          *(bcindJ+k)      = n2-2;
                          k++;                      }   
                         

  VecCreateSeqWithArray(MPI_COMM_SELF,4*localsizex,cacheScalar,&ybc);
  ISCreateGeneralWithIJ(MPI_COMM_SELF,*x,xvec,nvec,n2,4*localsizex, bcindI, bcindJ,ISfrom, ISto);
  LargeVecScatterCreate(xvec,ISfrom,ybc,ISto ,ctxt,nvec); 
  ISArrayDestroy(ISfrom,nvec);
  ISArrayDestroy(ISto,nvec);
  LargeVecScatterBeginEnd(xvec,ybc,INSERT_VALUES,SCATTER_FORWARD,ctxt,nvec);
  VecScatterArrayDestroy(ctxt,nvec);


  VecGetArray2d(ybc,4,localsizex,0,0,&ybcpt2);
  

  for(j=2;j<n2-2;j++){
     for(i=0;i<localsizex;i++){
      lvecptx[i][j] = C0*lvecpty[i][j]
                     +C1*lvecpty[i][j-1] +C1*lvecpty[i][j+1]
                     +C2*lvecpty[i][j-2] +C2*lvecpty[i][j+2];               
     }     
  }

  for(i=0;i<localsizex;i++){

      lvecptx[i][0] =  C0*lvecpty[i][0]
                       +C1*ybcpt2[1][i] +C1*lvecpty[i][1]
                       +C2*ybcpt2[0][i] +C2*lvecpty[i][2];  
      lvecptx[i][1] =  C0*lvecpty[i][1]
                       +C1*lvecpty[i][0] +C1*lvecpty[i][2]
                       +C2*ybcpt2[1][i] +C2*lvecpty[i][3]; 

      lvecptx[i][n2-2] = C0*lvecpty[i][n2-2]
                        +C1*lvecpty[i][n2-3] +C1*lvecpty[i][n2-1]
                        +C2*lvecpty[i][n2-4] +C2*ybcpt2[2][i];  

      lvecptx[i][n2-1] = C0*lvecpty[i][n2-1]
                        +C1*lvecpty[i][n2-2] +C1*ybcpt2[2][i]
                        +C2*lvecpty[i][n2-3] +C2*ybcpt2[3][i];  

  }

  VecRestoreArray2d(ybc,4,localsizex,0,0,&ybcpt2);
  VecDestroy(ybc);


/////////////////////////////////////////////////////////////////////////////////

  VecRestoreArray2d(*x,localsizex,n2,0,0,&lvecptx);
  VecRestoreArray2d(*y,localsizex,n2,0,0,&lvecpty);

  VecArrayDestroy(xvec,nvec); 




return 0;
}


//////////////////////////////////////////////////////////////////////////////////////

/*
Inverse Standard Map

Tzu-Chen Liang  11-6-2006 
*/
int InverseStandardMap(PetscScalar *x,PetscScalar *y, PetscScalar c){
 
  //PetscScalar c;
  //c = 0.3;

  *x = fmod(*x-*y,1);
  *y = fmod(*y-c*sin(2*M_PI*(*x)),1);

  if(*x<0){*x=*x+1;}
  if(*y<0){*y=*y+1;}

return 0;
}
//////////////////////////////////////////////////////////////////////////////////////
/*
Standard Map

Tzu-Chen Liang  11-8-2006
*/
int StandardMap(PetscScalar *x,PetscScalar *y, PetscScalar c){
 
  //PetscScalar c;
  //c = 0.3;

  *y = fmod(*y+c*sin(2*M_PI*(*x)),1);
  *x = fmod(*x+*y,1);

  if(*x<0){*x=*x+1;}
  if(*y<0){*y=*y+1;}
 

return 0;
}

//////////////////////////////////////////////////////////////////////////////////////

/*
Inverse InverseModifiedArnoldsCat Map

Tzu-Chen Liang  1-22-2007 
*/
int InverseModifiedArnoldsCatMap(PetscScalar *x,PetscScalar *y){
 
  PetscScalar c;

  c = (1e-3)/2/M_PI;

  *x = fmod(*x-*y,1);
  *y = fmod(*y-*x -c*sin(2*M_PI*(*x)),1);

  if(*x<0){*x=*x+1;}
  if(*y<0){*y=*y+1;}

return 0;
}


////////////////////////////////////////////////////////////////////////////
/*
Map a point (x,y)  to (1-x,1-y);

Tzu-Chen Liang  11-6-2006 
*/
int SkewSymmetricPoint(PetscInt *i, PetscInt *j, PetscInt n){
   *i = n - *i - 1;
   *j = n - *j - 1;
return 0;
}
////////////////////////////////////////////////////////////////////////////
int SkewSymmetricScatter(Vec *x,PetscScalar *cacheScalar ,PetscInt *cacheInt, PetscInt n2,PetscInt Istart,PetscInt localsizex , VecScatter *ctx){

PetscInt i,k; 
IS       isbc;
Vec      bcvec;

VecCreateSeqWithArray(MPI_COMM_SELF,localsizex*4,cacheScalar,&bcvec);

k = 0; 
for(i=0;i<localsizex;i++){*(cacheInt+k)= n2*(n2*2-Istart)-1-i*n2;k++;} 
for(i=0;i<localsizex;i++){*(cacheInt+k)= n2*(n2*2-Istart)-2-i*n2;k++;}
for(i=0;i<localsizex;i++){*(cacheInt+k)= n2*(n2*2-Istart)-n2+1-i*n2;k++;}
for(i=0;i<localsizex;i++){*(cacheInt+k)= n2*(n2*2-Istart)-n2-i*n2;k++;}

ISCreateGeneralWithArray(MPI_COMM_WORLD,4*localsizex,cacheInt,&isbc);
VecScatterCreate(*x,isbc,bcvec,PETSC_NULL,ctx);
ISDestroy(isbc);

return 0;
}





