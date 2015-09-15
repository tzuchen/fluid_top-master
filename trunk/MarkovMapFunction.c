#include "petscksp.h"
#include "petscda.h"
#include "slepceps.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "MeshStruc.h"
#include "PetscStreamline.h"
#include "PetscReceive.h"
#include "mpi.h"

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
int MarkovMap(Vec sol, PetscInt ny, PetscInt nz, CoordInfo *cdinfo,Mat *MapA, Mat6 *Dmat, Vec6 *vecf, PetscScalar **carray, PetscInt *nnzMapA)
{
    
    Vec            ix,iy,iz,jy,jz,tempx;
    Vec            j01y,j01z,j10y,j10z;
    Vec            q1,q2,ybar,zbar;
    Vec            d21,d2;

    Mat            A,Dyf,Dzf,Dy10f,Dz10f,Dy01f,Dz01f;

    PetscInt       Istart,Iend,Lsize,solsize;
    PetscScalar    yf,zf;
    PetscScalar    p11,p12,p21,p22;
    PetscScalar    detp, ip11,ip12,ip21,ip22;
    PetscScalar    dyz = 0.04/2,rowsum;
    PetscScalar    nzArray[140];
  
    PetscInt       i,j;
    PetscInt       Iind[1];
    PetscInt       nzcount, nzInd[140];
    
 
    PetscErrorCode ierr;
    PetscScalar    gsizey,gsizez;
    PetscInt       np, np2; 
  

    PetscInt n = 500;
  
    np   = ny*nz;
    np2  = np*np;

    PetscScalar x0[ny][nz],y0[ny][nz],z0[ny][nz];   
    PetscInt    dInd[np];

  
    for(i=0;i<np;i++){dInd[i]=i;}
    solsize   = cdinfo->nofgrid->x + cdinfo->nofgrid->y + cdinfo->nofgrid->z + 24;
    gsizey    = (PetscScalar)(cdinfo->r->y-cdinfo->l->y)/ny;
    gsizez    = (PetscScalar)(cdinfo->r->z-cdinfo->l->z)/nz;

    for(i=0;i<ny;i++){
         for (j=0;j<nz;j++){
                x0[i][j] = 0;
                y0[i][j] = (0.5+i)*gsizey;
                z0[i][j] = (0.5+j)*gsizez;      
         }
    }

    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)y0,&jy);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)z0,&jz);CHKERRQ(ierr);

    ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,np,&ix);CHKERRQ(ierr);
    ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,np,&iy);CHKERRQ(ierr);
    ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,np,&iz);CHKERRQ(ierr);

    ierr  = VecGetOwnershipRange(iy,&Istart,&Iend);CHKERRQ(ierr); 
    Lsize = Iend - Istart;

    PetscScalar zbarray[Lsize*140];
    PetscScalar ztarray[Lsize*140];
    PetscInt  Ind[Lsize];  
    PetscScalar ybarray[Lsize*140];
    PetscScalar ytarray[Lsize*140];
  //a fake one due to memory free problem!
    PetscScalar yytarray[Lsize*140];

    PetscInt nnzcount = 0;
    PetscScalar *darray;

  
    for(i=0;i<Lsize;i++){Ind[i]=Istart+i;}
    ierr  =  VecSetValues(ix,Lsize,Ind,(PetscScalar*)x0+Istart,INSERT_VALUES);CHKERRQ(ierr);
    ierr  =  VecSetValues(iy,Lsize,Ind,(PetscScalar*)y0+Istart,INSERT_VALUES);CHKERRQ(ierr);
    ierr  =  VecSetValues(iz,Lsize,Ind,(PetscScalar*)z0+Istart,INSERT_VALUES);CHKERRQ(ierr);
   
    ierr  =  VecAssemblyBegin(ix);CHKERRQ(ierr);
    ierr  =  VecAssemblyEnd(ix);CHKERRQ(ierr);
    ierr  =  VecAssemblyBegin(iy);CHKERRQ(ierr);
    ierr  =  VecAssemblyEnd(iy);CHKERRQ(ierr);
    ierr  =  VecAssemblyBegin(iz);CHKERRQ(ierr); 
    ierr  =  VecAssemblyEnd(iz);CHKERRQ(ierr);

    
    ierr  = VecDuplicate(iy,&j01y);CHKERRQ(ierr); 
    ierr  = VecDuplicate(iz,&j01z);CHKERRQ(ierr); 
    ierr  = VecDuplicate(iy,&j10y);CHKERRQ(ierr); 
    ierr  = VecDuplicate(iz,&j10z);CHKERRQ(ierr); 


    ierr  = VecDuplicate(jy,&d21);CHKERRQ(ierr); 
    ierr  = VecDuplicate(jz,&d2);CHKERRQ(ierr); 
    ierr  = VecDuplicate(jy,&q1);CHKERRQ(ierr); 
    ierr  = VecDuplicate(jz,&q2);CHKERRQ(ierr); 
    ierr  = VecDuplicate(jy,&ybar);CHKERRQ(ierr); 
    ierr  = VecDuplicate(jz,&zbar);CHKERRQ(ierr); 

    ierr  = VecCopy(iy,j01y);CHKERRQ(ierr); 
    ierr  = VecCopy(iz,j01z);CHKERRQ(ierr);
    ierr  = VecCopy(iy,j10y);CHKERRQ(ierr);
    ierr  = VecCopy(iz,j10z);CHKERRQ(ierr);

    ierr  = VecShift(j10y,dyz);CHKERRQ(ierr);
    ierr  = VecShift(j01z,dyz);CHKERRQ(ierr);

    ierr  = VecDuplicate(ix,&tempx);CHKERRQ(ierr);
    StreamlineMapr(&tempx, &iy,   &iz,   sol, cdinfo, &Dyf, &Dzf,32, n ,10);
    ierr  = VecDuplicate(ix,&tempx);CHKERRQ(ierr);
    StreamlineMapr(&tempx, &j10y, &j10z, sol, cdinfo, &Dy10f, &Dz10f,32, n ,10);
    ierr  = VecDuplicate(ix,&tempx);CHKERRQ(ierr);
    StreamlineMapr(&tempx, &j01y, &j01z, sol, cdinfo, &Dy01f, &Dz01f,32, n ,10);
 
    MatCreateMPIAIJ(MPI_COMM_WORLD,Lsize,np,PETSC_DETERMINE,PETSC_DETERMINE,140,PETSC_NULL,140,PETSC_NULL,&A);CHKERRQ(ierr);

    for(i=Istart;i<Iend;i++){
        Iind[0] = i;
   
       ierr  = VecGetValues(iy,1,Iind,&yf);CHKERRQ(ierr);
       ierr  = VecGetValues(iz,1,Iind,&zf);CHKERRQ(ierr);       
       ierr  = VecGetValues(j10y,1,Iind,&p11);CHKERRQ(ierr);
       ierr  = VecGetValues(j10z,1,Iind,&p21);CHKERRQ(ierr);
       ierr  = VecGetValues(j01y,1,Iind,&p12);CHKERRQ(ierr);
       ierr  = VecGetValues(j01z,1,Iind,&p22);CHKERRQ(ierr);  

       p11   = 1;//(p11 - yf)/dyz;
       p21   = 0;//(p21 - zf)/dyz;
       p12   = 0;//(p12 - yf)/dyz;
       p22   = 1;//(p22 - zf)/dyz;     

       detp  =  p11*p22-p21*p12;
       ip11  =  p22/detp;
       ip12  = -p12/detp;
       ip21  = -p21/detp;
       ip22  =  p11/detp;
  
       ierr  = VecCopy(jy,q1);CHKERRQ(ierr);    
       ierr  = VecCopy(jz,q2);CHKERRQ(ierr);    
       ierr  = VecShift(q1,-yf);CHKERRQ(ierr);
       ierr  = VecShift(q2,-zf);CHKERRQ(ierr);

       ierr  = VecCopy(q2,ybar);CHKERRQ(ierr);
       ierr  = VecCopy(q2,zbar);CHKERRQ(ierr);
       ierr  = VecAXPBY(ybar,ip11,ip12,q1);
       ierr  = VecAXPBY(zbar,ip21,ip22,q1);

       ierr  = VecPointwiseMult(d21,ybar,ybar);CHKERRQ(ierr);
       ierr  = VecPointwiseMult(d2 ,zbar,zbar);CHKERRQ(ierr);
       ierr  = VecAXPY(d2,1,d21);CHKERRQ(ierr);
     
       ierr =  VecGetArray(d2,&darray);CHKERRQ(ierr);
       nzcount = 0;
       rowsum  = 0;
       for(j=0;j<np;j++) { 
           if(darray[j]<25*gsizey*gsizez && nzcount<140){

                  nzInd[nzcount]   = j;
                  nzArray[nzcount] = exp(-darray[j]/gsizey/gsizez);
                  rowsum           = rowsum + nzArray[nzcount];
                  nzcount          = nzcount+1;    
                  
                }else{
                  darray[j] = 0;                   
                }
        }
        VecRestoreArray(d2,&darray);    
        VecGetValues(ybar,nzcount,nzInd,ybarray+nnzcount);
        VecGetValues(zbar,nzcount,nzInd,zbarray+nnzcount);

        for(j=nnzcount;j<nnzcount+nzcount;j++){
             ytarray[j] = ip11*ybarray[j] + ip21*zbarray[j];
             ztarray[j] = ip12*ybarray[j] + ip22*zbarray[j];
        } 
        nnzcount = nnzcount + nzcount;

        for (j=0;j<nzcount;j++){nzArray[j] = nzArray[j]/rowsum;}   
  
        // ierr  = VecSetValues(d2,np,dInd,darray,INSERT_VALUES);CHKERRQ(ierr);  
        ierr  = MatSetValues(A,1,Iind,nzcount,nzInd,nzArray, INSERT_VALUES);CHKERRQ(ierr);  
   }


   ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
   ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);

   PetscScalar c1[nnzcount],c2[nnzcount],c3[nnzcount],c4[nnzcount],c5[nnzcount],c6[nnzcount];


   for (j=0;j<nnzcount;j++){
          c3[j] =  0*ytarray[j]*ybarray[j];
          c4[j] =  0*ztarray[j]*ybarray[j];   
          c5[j] =  0*ytarray[j]*zbarray[j];
          c6[j] =  0*ztarray[j]*zbarray[j];
          c1[j] = -c3[j]-c5[j]-ytarray[j];
          c2[j] = -c4[j]-c6[j]-ztarray[j];
   }


   
    carray[0]   = c1;
    carray[1]   = c2;
    carray[2]   = c3;
    carray[3]   = c4;
    carray[4]   = c5;
    carray[5]   = c6;
  
   *nnzMapA      = nnzcount;
   *MapA         = A;

    Dmat->Mat1 = Dyf;
    Dmat->Mat2 = Dzf;
    Dmat->Mat3 = Dy10f;
    Dmat->Mat4 = Dz10f;
    Dmat->Mat5 = Dy01f;
    Dmat->Mat6 = Dz01f;

    vecf->vec1 = iy;
    vecf->vec2 = iz;
    vecf->vec3 = j10y;
    vecf->vec4 = j10z;
    vecf->vec5 = j01y;
    vecf->vec6 = j01z;
  

    free(zbarray);
    free(ztarray);
    free(ybarray);
    free(ytarray);

 // free(yytarray);


    free(Ind);
 // free(dInd);

    ierr = VecDestroy(q1);CHKERRQ(ierr);
    ierr = VecDestroy(q2);CHKERRQ(ierr);
    ierr = VecDestroy(d21);CHKERRQ(ierr);
    ierr = VecDestroy(d2);CHKERRQ(ierr);
    ierr = VecDestroy(ix);CHKERRQ(ierr);
    //ierr = VecDestroy(ybar);CHKERRQ(ierr);
    //ierr = VecDestroy(zbar);CHKERRQ(ierr);
 
return 0;

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
int dMarkovMap(Vec *dAvec, Vec6 cvec, Mat6 Dmat)
{

  Vec            dA;
  PetscInt       m,n;
  PetscErrorCode ierr;

  MatGetSize(Dmat.Mat1,&m,&n);

  ierr  = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,n,&dA);CHKERRQ(ierr);

  ierr  = MatMultTranspose(Dmat.Mat1,cvec.vec1,dA);CHKERRQ(ierr);
  ierr  = MatMultTransposeAdd(Dmat.Mat2,cvec.vec2,dA,dA);CHKERRQ(ierr);
  ierr  = MatMultTransposeAdd(Dmat.Mat3,cvec.vec3,dA,dA);CHKERRQ(ierr);
  ierr  = MatMultTransposeAdd(Dmat.Mat4,cvec.vec4,dA,dA);CHKERRQ(ierr);
  ierr  = MatMultTransposeAdd(Dmat.Mat5,cvec.vec5,dA,dA);CHKERRQ(ierr);
  ierr  = MatMultTransposeAdd(Dmat.Mat6,cvec.vec6,dA,dA);CHKERRQ(ierr);

  *dAvec = dA;
 return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// 9pt method, more memory required
// 
int MarkovMap2(Vec sol, PetscInt ny, PetscInt nz, CoordInfo *cdinfo,Mat *MapA)
{

    Vec         temp,x0vec,y0vec,z0vec;
    Mat         Aseq,A;
   
     
    PetscInt    yStart,yEnd,ySize;
    PetscInt    ind,np,n,i,j,*col,p,q,k,l;
    PetscInt    *iy,*iz,*yshift,*zshift;
    PetscScalar *x0, *y0, *z0,*val;
    PetscScalar gsizey,gsizez,hgsizey,hgsizez,valsum;
    PetscScalar *ye,*ze,py,pz;
    PetscErrorCode ierr;


    gsizey = (cdinfo->r->y-cdinfo->l->y)/(1.0*ny);
    gsizez = (cdinfo->r->z-cdinfo->l->z)/(1.0*nz);
    hgsizey = gsizey/2;
    hgsizez = gsizez/2;
    n      = 500;

     
    VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,ny,&temp);
    VecGetOwnershipRange(temp,&yStart,&yEnd);    
    ySize = yEnd-yStart;

    np = (2*ySize+1)*(2*nz+1);
    x0 = malloc(np*sizeof(PetscScalar));
    y0 = malloc(np*sizeof(PetscScalar));
    z0 = malloc(np*sizeof(PetscScalar));

    ye = malloc(9*sizeof(PetscScalar));
    ze = malloc(9*sizeof(PetscScalar));
    val = malloc(9*sizeof(PetscScalar));

    yshift = malloc(9*sizeof(PetscInt));
    zshift = malloc(9*sizeof(PetscInt));
    iy = malloc(9*sizeof(PetscInt));
    iz = malloc(9*sizeof(PetscInt));
    col = malloc(9*sizeof(PetscInt));

    ind = 0;
    for(i=-1;i<2;i++){
       for(j=-1;j<2;j++){
           *(yshift+ind) = i;
           *(zshift+ind) = j;
           ind = ind+1;
          }
     }


    ind = 0;
    for(i=2*yStart;i<=2*yEnd;i++){
         for (j=0;j<=2*nz;j++){
                *(x0+ind) = 0;
                *(y0+ind) = i*hgsizey;
                *(z0+ind) = j*hgsizez;      
                ind = ind+1;
         }
    }
    // use seq vector
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)x0,&x0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)y0,&y0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)z0,&z0vec);CHKERRQ(ierr);
   
    StreamlineMaprr(&x0vec, &y0vec,   &z0vec,   sol, cdinfo,32, n ,10);
  

    MatCreateSeqAIJ(PETSC_COMM_SELF,ySize*nz,ny*nz,9,PETSC_NULL,&Aseq);
  
    ind = 0;
        for(i=0;i<ySize;i++){
          for (j=0;j<nz;j++){
             
              p = 2*i+1;
              q = 2*j+1;
 
              for(k=0;k<9;k++){
                p   = 2*i+1+*(yshift+k);
                q   = 2*j+1+*(zshift+k);
                *(ye+k) = *(y0+p*(2*nz+1)+q);
                *(ze+k) = *(z0+p*(2*nz+1)+q); 
                *(iy+k) = (PetscInt)floor(*(ye+k)/gsizey);
                *(iz+k) = (PetscInt)floor(*(ze+k)/gsizez);   

                 if(*(iy+k)>=ny){*(iy+k)=ny-1;}
                 if(*(iz+k)>=ny){*(iz+k)=nz-1;}
                 if(*(iy+k)<0){*(iy+k)=0;}
                 if(*(iz+k)<0){*(iz+k)=0;}

                *(col+k) = *(iy+k)*nz+*(iz+k);
                   for(l=0;l<k;l++){
                      if(*(col+l)==*(col+k)) {*(col+k) =-1;} 
                   }
              }
  
               valsum = 0; 
               for(k=0;k<9;k++){
                   if(*(col+k)>0){
                      py = fabs(*(iy+k)*gsizey+hgsizey-*(ye+4));
                      pz = fabs(*(iz+k)*gsizez+hgsizez-*(ze+4));
        
                      *(val+k)=fabs(gsizey-py)*fabs(gsizez-pz); 
                       valsum = valsum+*(val+k);
                    }else{
                      *(val+k) = 0;
                    }
               }
               for(k=0;k<9;k++){
                  *(val+k)= *(val+k)/valsum;
               }
            
       
              MatSetValues(Aseq,1,&ind,9,col,val,ADD_VALUES);
           
             ind = ind+1;
          }
        }

   MatAssemblyBegin(Aseq,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(Aseq,MAT_FINAL_ASSEMBLY); 

        MatMerge(MPI_COMM_WORLD,Aseq,ySize*nz,MAT_INITIAL_MATRIX,&A);
        *MapA = A; 

////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////


        ierr = VecDestroy(x0vec);CHKERRQ(ierr);
        ierr = VecDestroy(y0vec);CHKERRQ(ierr);
        ierr = VecDestroy(z0vec);CHKERRQ(ierr);
        free(ye);free(ze);free(val); 
        free(yshift); free(zshift);
        free(iy);free(iz);free(col);
        
  
return 0;

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// Now try 9pt weighted method, hm..
int MarkovMap3(Vec sol, PetscInt ny, PetscInt nz, CoordInfo *cdinfo,Mat *MapA)
{
 
    Vec         temp,x0vec,y0vec,z0vec;
    Mat         Aseq,A;
   
     
    PetscInt    yStart,yEnd,ySize;
    PetscInt    ind,np,n,i,j,k,p;
    PetscInt    *iy,*iz,*col;
    PetscScalar *x0, *y0, *z0;
    PetscScalar *xe,*ye,*ze,*val,*yshift,*zshift,*dyshift,*dzshift;
    PetscScalar gsizey,gsizez,hgsizey,hgsizez;
    PetscScalar py,pz,valsum;
    PetscErrorCode ierr;


    gsizey = (cdinfo->r->y-cdinfo->l->y)/(1.0*ny);
    gsizez = (cdinfo->r->z-cdinfo->l->z)/(1.0*nz);
    n      = 500;

    hgsizey = gsizey/2;
    hgsizez = gsizez/2;
 
  
     
    VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,ny,&temp);
    VecGetOwnershipRange(temp,&yStart,&yEnd);    
    ySize = yEnd-yStart;

    

    np = ySize*nz;
    x0 = malloc(np*sizeof(PetscScalar));
    y0 = malloc(np*sizeof(PetscScalar));
    z0 = malloc(np*sizeof(PetscScalar));
 
    xe = malloc(9*sizeof(PetscScalar));
    ye = malloc(9*sizeof(PetscScalar));
    ze = malloc(9*sizeof(PetscScalar));
    yshift = malloc(9*sizeof(PetscScalar));
    zshift = malloc(9*sizeof(PetscScalar));
    dyshift = malloc(9*sizeof(PetscScalar));
    dzshift = malloc(9*sizeof(PetscScalar));

    iy = malloc(9*sizeof(PetscInt));
    iz = malloc(9*sizeof(PetscInt));
   col = malloc(9*sizeof(PetscInt));
   val = malloc(9*sizeof(PetscScalar));

    for(i=0;i<9;i++){*(val+i)=1.0/9;}

    ind = 0;
    for(i=yStart;i<yEnd;i++){
         for (j=0;j<nz;j++){
                *(x0+ind) = 0;
                *(y0+ind) = i*gsizey+hgsizey;
                *(z0+ind) = j*gsizez+hgsizez;
      
                ind = ind+1;
         }
    }

    ind = 0;
    for(i=-1;i<2;i++){
         for (j=-1;j<2;j++){
          *(yshift+ind) = i*(hgsizey-0*1e-5);
          *(zshift+ind) = j*(hgsizez-0*1e-5);
          *(dyshift+ind) = -i*1e-5;
          *(dzshift+ind) = -j*1e-5;
          ind=ind+1;
         }
    }
  

    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,9,(PetscScalar*)xe,&x0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,9,(PetscScalar*)ye,&y0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,9,(PetscScalar*)ze,&z0vec);CHKERRQ(ierr);
   
 

    MatCreateSeqAIJ(PETSC_COMM_SELF,ySize*nz,ny*nz,9,PETSC_NULL,&Aseq);
  
    ind = 0;
        for(i=0;i<ySize;i++){
          for (j=0;j<nz;j++){
   
              for(k=0;k<9;k++){
                    *(xe+k) = 0;
                    *(ye+k) = *(y0+ind) + *(yshift+k);
                    *(ze+k) = *(z0+ind) + *(zshift+k);
                    
              }

        // StreamlineMapr(&x0vec, &y0vec,   &z0vec,   sol, cdinfo, &Dyf, &Dzf,32, n ,10);
         StreamlineMaprr(&x0vec, &y0vec,   &z0vec,   sol, cdinfo,32, n ,10);

              for(k=0;k<9;k++){
                    *(ye+k) = *(ye+k)  + *(dyshift+k);
                    *(ze+k) = *(ze+k)  + *(dzshift+k);

                    *(iy+k) = (PetscInt)floor(*(ye+k)/gsizey); 
                    *(iz+k) = (PetscInt)floor(*(ze+k)/gsizez); 
                    *(col+k) = *(iy+k)*nz+*(iz+k);
 
                    for(p=0;p<k;p++){
                      if(*(col+p)==*(col+k)) {*(col+k) =-1;} 
                   }

              }



        


               valsum = 0; 
               for(k=0;k<9;k++){
                   if(*(col+k)>0){
                      py = fabs(*(iy+k)*gsizey+hgsizey-*(ye+4));
                      pz = fabs(*(iz+k)*gsizez+hgsizez-*(ze+4));
        
                      *(val+k)=fabs(gsizey-py)*fabs(gsizez-pz); 
                       valsum = valsum+*(val+k);
                    }else{
                      *(val+k) = 0;
                    }
               }
               for(k=0;k<9;k++){
                  *(val+k)= *(val+k)/valsum;
               }

               MatSetValues(Aseq,1,&ind,9,col,val,ADD_VALUES);
               ind=ind+1;

          }
        }
 
        MatMerge(MPI_COMM_WORLD,Aseq,ySize*nz,MAT_INITIAL_MATRIX,&A);
        *MapA = A; 

        ierr = VecDestroy(x0vec);CHKERRQ(ierr);
        ierr = VecDestroy(y0vec);CHKERRQ(ierr);
        ierr = VecDestroy(z0vec);CHKERRQ(ierr);
      //  ierr = MatDestroy(Aseq);CHKERRQ(ierr);




        free(x0); free(y0); free(z0);
        free(xe); free(ye); free(ze);
        free(yshift); free(zshift);
        free(dyshift); free(dzshift);
        free(iy); free(iz);
        free(col); free(val);
      
  



return 0;

}
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
// Sample method, more memory required
int MarkovMap4(Vec sol, PetscInt ny, PetscInt nz, CoordInfo *cdinfo,Mat *MapA)
{

  Vec            temp;
  Vec            x0vec, y0vec, z0vec;
  Mat            A,Aseq;
  PetscErrorCode ierr;

  PetscScalar       gsizey,gsizez,hgsizey,hgsizez;
  PetscInt       den,n,np,ind,i,j,p,q;
  PetscInt       ySize, yStart, yEnd;
  PetscScalar    *x0,*y0,*z0,*val;
  PetscInt       *iy,*iz,*col,*row;
  
 

  ierr  = VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,ny,&temp);CHKERRQ(ierr);
  ierr  = VecGetOwnershipRange(temp,&yStart,&yEnd);CHKERRQ(ierr);    
  ySize = yEnd-yStart;


  n  =500;
  den = 5;

  gsizey = (cdinfo->r->y-cdinfo->l->y)/(1.0*ny);
  gsizez = (cdinfo->r->z-cdinfo->l->z)/(1.0*nz);
  hgsizey = gsizey/(den-1);
  hgsizez = gsizez/(den-1);

  np  = (den*ySize)*(den*nz);
  x0 = malloc(np*sizeof(PetscScalar));
  y0 = malloc(np*sizeof(PetscScalar));
  z0 = malloc(np*sizeof(PetscScalar));
  
  iy  = malloc(np*sizeof(PetscInt));
  iz  = malloc(np*sizeof(PetscInt));
  col = malloc(np*sizeof(PetscInt));
  val = malloc(np*sizeof(PetscInt));
  row = malloc(np*sizeof(PetscInt));
 


  ind = 0;
  for(i=yStart;i<yEnd;i++){
      for (j=0;j<nz;j++){
           for(p=0;p<den;p++){           
                for(q=0;q<den;q++){
                  *(x0+ind) = 0;
                  //*(y0+ind) = i*gsizey+(p+0.5)*hgsizey;
                  //*(z0+ind) = j*gsizez+(q+0.5)*hgsizez;    
                   *(y0+ind) = i*gsizey+p*hgsizey;
                   *(z0+ind) = j*gsizez+q*hgsizez; 
                  *(row+ind) = (i-yStart)*nz+j;            
                   ind = ind+1;
               }
           }
       }
   }

    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)x0,&x0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)y0,&y0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)z0,&z0vec);CHKERRQ(ierr);   
    StreamlineMaprr(&x0vec, &y0vec,   &z0vec,   sol, cdinfo,32, n ,10);

    for(ind=0;ind<np;ind++){
            *(iy+ind)  = (PetscInt)floor(*(y0+ind)/gsizey);
            *(iz+ind)  = (PetscInt)floor(*(z0+ind)/gsizez);  
                
                 if(*(iy+ind)>=ny){*(iy+ind)=ny-1;}
                 if(*(iz+ind)>=ny){*(iz+ind)=nz-1;}
                 if(*(iy+ind)<0){*(iy+ind)=0;}
                 if(*(iz+ind)<0){*(iz+ind)=0;}

 
            *(col+ind) = *(iy+ind)*nz+*(iz+ind);
            *(val+ind) = (1.0/den)/den;

    }

 
   ierr  = MatCreateSeqAIJ(PETSC_COMM_SELF,ySize*nz,ny*nz,9,PETSC_NULL,&Aseq);CHKERRQ(ierr);

   for(i=0;i<ySize*nz;i++){
   ierr  = MatSetValues(Aseq,1,&i,den*den,(col+i*den*den),val,ADD_VALUES);CHKERRQ(ierr);
   }
 
   MatAssemblyBegin(Aseq,MAT_FINAL_ASSEMBLY);
   MatAssemblyEnd(Aseq,MAT_FINAL_ASSEMBLY);



   ierr  = MatMerge(MPI_COMM_WORLD,Aseq,ySize*nz,MAT_INITIAL_MATRIX,&A);CHKERRQ(ierr);
   *MapA = A; 
        
  
return 0;

}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*
9pt weighted method, but this function returns dA/dv also

Tzu-Chen Liang 3-5-2007
*/
int MarkovMap2d(Vec sol, PetscInt ny, PetscInt nz, CoordInfo *cdinfo,Mat *MapA, Mat *dAdy, Mat *dAdz)
{

    Vec         temp,x0vec,y0vec,z0vec;
    Mat         Aseq,A;
    Mat         DAdy,DAdz,DAdyseq,DAdzseq;  
     
    PetscInt    yStart,yEnd,ySize;
    PetscInt    ind,np,n,i,j,*col,p,q,k,l;
    PetscInt    *iy,*iz,*yshift,*zshift;
    PetscScalar *x0, *y0, *z0,*val;
    PetscScalar *dAdyval, *dAdzval;
    PetscScalar gsizey,gsizez,hgsizey,hgsizez,valsum;
    PetscScalar *ye,*ze,py,pz,wy,wz;
    PetscErrorCode ierr;
    PetscInt     solmaxp;
    PetscScalar  solmax;
 

    gsizey = (cdinfo->r->y-cdinfo->l->y)/(1.0*ny);
    gsizez = (cdinfo->r->z-cdinfo->l->z)/(1.0*nz);



    hgsizey = gsizey/2.0;
    hgsizez = gsizez/2.0;
    n      = 2000;


    ierr = VecMax(sol,&solmaxp,&solmax);CHKERRQ(ierr);
     
    ierr = VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,ny,&temp);CHKERRQ(ierr);
    ierr = VecGetOwnershipRange(temp,&yStart,&yEnd); CHKERRQ(ierr);   
    ySize = yEnd-yStart;

    np      = (2*ySize+1)*(2*nz+1);
    x0      = malloc(np*sizeof(PetscScalar));
    y0      = malloc(np*sizeof(PetscScalar));
    z0      = malloc(np*sizeof(PetscScalar));

    ye      = malloc(9*sizeof(PetscScalar));
    ze      = malloc(9*sizeof(PetscScalar));
    val     = malloc(9*sizeof(PetscScalar));
    dAdyval = malloc(9*sizeof(PetscScalar));
    dAdzval = malloc(9*sizeof(PetscScalar));

    yshift  = malloc(9*sizeof(PetscInt));
    zshift  = malloc(9*sizeof(PetscInt));
    iy      = malloc(9*sizeof(PetscInt));
    iz      = malloc(9*sizeof(PetscInt));
    col     = malloc(9*sizeof(PetscInt));

    ind = 0;
    for(i=-1;i<2;i++){
       for(j=-1;j<2;j++){
           *(yshift+ind) = i;
           *(zshift+ind) = j;
           ind = ind+1;
          }
     }

    ind = 0;
    for(i=2*yStart;i<=2*yEnd;i++){
         for (j=0;j<=2*nz;j++){
                *(x0+ind) = 0;
                *(y0+ind) = i*hgsizey;
                *(z0+ind) = j*hgsizez;      
                ind = ind+1;
         }
    }
    // use seq vector (they have one row overlap)
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)x0,&x0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)y0,&y0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)z0,&z0vec);CHKERRQ(ierr);
    
    //dt = 0.01/solmax
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: solmax= %f \n",solmax);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"PETSC: dt = %f  \n",0.01/solmax);CHKERRQ(ierr); 

    StreamlineMaprr(&x0vec, &y0vec,   &z0vec,   sol, cdinfo,4*0.01/solmax,  n ,10);
    VecSave(y0vec,"yefiles");
    VecSave(z0vec,"zefiles");

    ierr  = MatCreateSeqAIJ(PETSC_COMM_SELF,ySize*nz,ny*nz,9,PETSC_NULL,&Aseq);CHKERRQ(ierr);
    ierr  = MatCreateSeqAIJ(PETSC_COMM_SELF,ySize*nz,ny*nz,9,PETSC_NULL,&DAdyseq);CHKERRQ(ierr);
    ierr  = MatCreateSeqAIJ(PETSC_COMM_SELF,ySize*nz,ny*nz,9,PETSC_NULL,&DAdzseq);CHKERRQ(ierr);
  
        ind = 0;
        for(i=0;i<ySize;i++){
              for (j=0;j<nz;j++){
             
                  p = 2*i+1;
                  q = 2*j+1;
 
                  for(k=0;k<9;k++){
                     p   = 2*i+1+*(yshift+k);
                     q   = 2*j+1+*(zshift+k);
                     *(ye+k) = *(y0+p*(2*nz+1)+q);
                     *(ze+k) = *(z0+p*(2*nz+1)+q); 
                     *(iy+k) = (PetscInt)floor(*(ye+k)/gsizey);
                     *(iz+k) = (PetscInt)floor(*(ze+k)/gsizez);   

                      //if(*(iy+k)>=ny){*(iy+k)=ny-1;}
                      //if(*(iz+k)>=nz){*(iz+k)=nz-1;}
                      //if(*(iy+k)<0){*(iy+k)=0;}
                      //if(*(iz+k)<0){*(iz+k)=0;}

                     if(cdinfo->period->y==0){
                      if(*(iy+k)>=ny){*(iy+k)=ny-1;}
                      if(*(iy+k)<0){*(iy+k)=0;}
                     }else{
                      if(*(iy+k)>=ny){*(iy+k)=*(iy+k)-ny;}
                      if(*(iy+k)<0){*(iy+k)=*(iy+k)+ny;}
                     }

                     if(cdinfo->period->z==0){
                         if(*(iz+k)>=nz){*(iz+k)=nz-1;}                      
                         if(*(iz+k)<0){*(iz+k)=0;}
                     }else{
                         if(*(iz+k)>=nz){*(iz+k)=*(iz+k)-nz;}                      
                         if(*(iz+k)<0){*(iz+k)=*(iz+k)+nz;}
                      }

 
                     *(col+k) = *(iy+k)*nz+*(iz+k);
                     for(l=0;l<k;l++){ if(*(col+l)==*(col+k)) {*(col+k) =-1;}}
                  }
  
                  valsum = 0; 
                  for(k=0;k<9;k++){
                      if(*(col+k)>0){
                         py = *(ye+4)-(*(iy+k)*gsizey+hgsizey);
                         pz = *(ze+4)-(*(iz+k)*gsizez+hgsizez);
                         
                         wy = gsizey-fabs(py);
                         if(wy<0){ wy=0;}
                         wz = gsizez-fabs(pz);
                         if(wz<0){ wz=0;}
                         *(val+k)=wy*wz;
                         if(py>0){*(dAdyval+k) =  wz;}
                             else{*(dAdyval+k) = -wz;} 
                         if(pz>0){*(dAdzval+k) =  wy;}
                             else{*(dAdzval+k) = -wy;} 
                            
                         //*(val+k)=fabs(gsizey-fabs(py))*fabs(gsizez-fabs(pz));
                         
                        // if(py>0){*(dAdyval+k) =  fabs(gsizez-fabs(pz));}
                        //     else{*(dAdyval+k) = -fabs(gsizez-fabs(pz));} 
                        // if(pz>0){*(dAdzval+k) =  fabs(gsizey-fabs(py));}
                        //    else{*(dAdzval+k) = -fabs(gsizey-fabs(py));} 
                         valsum = valsum+*(val+k);
                      }else{
                         *(val+k) = 0;
                      }
                  }
                  for(k=0;k<9;k++){*(val+k)= *(val+k)/valsum;if(isnan(*(val+k))){*(val+k)=1;}  }
            
                  ierr = MatSetValues(Aseq,1,&ind,9,col,val,ADD_VALUES);CHKERRQ(ierr);
                  ierr = MatSetValues(DAdyseq,1,&ind,9,col,dAdyval,ADD_VALUES);CHKERRQ(ierr);
                  ierr = MatSetValues(DAdzseq,1,&ind,9,col,dAdzval,ADD_VALUES);CHKERRQ(ierr);
                  ind = ind+1;
             }
        }
      

        ierr = MatAssemblyBegin(Aseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Aseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
        ierr = MatAssemblyBegin(DAdyseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(DAdyseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyBegin(DAdzseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(DAdzseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

        ierr =  MatMerge(MPI_COMM_WORLD,Aseq,ySize*nz,MAT_INITIAL_MATRIX,&A);CHKERRQ(ierr);
        ierr =  MatMerge(MPI_COMM_WORLD,DAdyseq,ySize*nz,MAT_INITIAL_MATRIX,&DAdy);CHKERRQ(ierr);
        ierr =  MatMerge(MPI_COMM_WORLD,DAdzseq,ySize*nz,MAT_INITIAL_MATRIX,&DAdz);CHKERRQ(ierr);
        
                     //MatSetValue(A,0,0,1,INSERT_VALUES);

        *MapA = A; 
        *dAdy = DAdy;
        *dAdz = DAdz;
   

////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////


        ierr = VecDestroy(x0vec);CHKERRQ(ierr);
        ierr = VecDestroy(y0vec);CHKERRQ(ierr);
        ierr = VecDestroy(z0vec);CHKERRQ(ierr);
        free(ye);free(ze);free(val); 
        free(yshift); free(zshift);
        free(iy);free(iz);free(col);
        
  
return 0;

}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/*
9pt weighted method, but this function returns dA/dv also
And we normalize A by the velocity in x dir

Tzu-Chen Liang 9-12-2007
*/
int MarkovMap2dv(Vec sol, PetscInt ny, PetscInt nz, CoordInfo *cdinfo,Mat *MapA, Mat *dAdy, Mat *dAdz)
{

    Vec         temp,x0vec,y0vec,z0vec;
    Mat         Aseq,A;
    Mat         DAdy,DAdz,DAdyseq,DAdzseq;  
     
    PetscInt    yStart,yEnd,ySize;
    PetscInt    ind,np,n,i,j,*col,p,q,k,l;
    PetscInt    *iy,*iz,*yshift,*zshift;
    PetscScalar *x0, *y0, *z0,*val;
    PetscScalar *dAdyval, *dAdzval;
    PetscScalar gsizey,gsizez,hgsizey,hgsizez,valsum;
    PetscScalar *ye,*ze,py,pz;
    PetscErrorCode ierr;
    
     Vec            vx0vec;
     PetscScalar    *vx0array;
     PetscInt       pind[24];
     PetscScalar    weight[24];
     PetscScalar3   vel;
     PetscScalar    *velx0;
     PetscScalar    vr;

    gsizey = (cdinfo->r->y-cdinfo->l->y)/(1.0*ny);
    gsizez = (cdinfo->r->z-cdinfo->l->z)/(1.0*nz);
    hgsizey = gsizey/2;
    hgsizez = gsizez/2;
    n      = 5500;


   

    
    //ierr = VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,ny,&vx0vec);
     
    ierr = VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,ny,&temp);
    ierr = VecGetOwnershipRange(temp,&yStart,&yEnd);    
    ySize = yEnd-yStart;

    np      = (2*ySize+1)*(2*nz+1);
    x0      = malloc(np*sizeof(PetscScalar));
    y0      = malloc(np*sizeof(PetscScalar));
    z0      = malloc(np*sizeof(PetscScalar));

    ye      = malloc(9*sizeof(PetscScalar));
    ze      = malloc(9*sizeof(PetscScalar));
    val     = malloc(9*sizeof(PetscScalar));
    dAdyval = malloc(9*sizeof(PetscScalar));
    dAdzval = malloc(9*sizeof(PetscScalar));

    yshift  = malloc(9*sizeof(PetscInt));
    zshift  = malloc(9*sizeof(PetscInt));
    iy      = malloc(9*sizeof(PetscInt));
    iz      = malloc(9*sizeof(PetscInt));
    col     = malloc(9*sizeof(PetscInt));

    ind = 0;
    for(i=-1;i<2;i++){
       for(j=-1;j<2;j++){
           *(yshift+ind) = i;
           *(zshift+ind) = j;
           ind = ind+1;
          }
     }

    ind = 0;
    for(i=2*yStart;i<=2*yEnd;i++){
         for (j=0;j<=2*nz;j++){
                *(x0+ind) = 0;
                *(y0+ind) = i*hgsizey;
                *(z0+ind) = j*hgsizez;      
                ind = ind+1;
         }
    }
    // use seq vector
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)x0,&x0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)y0,&y0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)z0,&z0vec);CHKERRQ(ierr);
   


    k = 0; 
    velx0      = malloc(ny*nz*sizeof(PetscScalar));
    for(i=0;i<ny;i++){
      for(j=0;j<nz;j++){       
      veval(0, (i-0.5)*gsizey, (j-0.5)*gsizez, &sol, cdinfo, &vel,pind, weight);  
     *(velx0+k) = vel.x;
     k++;
      }
    }
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,ny*nz,(PetscScalar*)velx0,&vx0vec);CHKERRQ(ierr);
    ierr = VecSave(vx0vec,"vx0vecfiles");CHKERRQ(ierr);


    StreamlineMaprr(&x0vec, &y0vec,   &z0vec,   sol, cdinfo,32,  n ,10);

    VecSave(y0vec,"yefiles");
    VecSave(z0vec,"zefiles");

    ierr  = MatCreateSeqAIJ(PETSC_COMM_SELF,ySize*nz,ny*nz,9,PETSC_NULL,&Aseq);CHKERRQ(ierr);
    ierr  = MatCreateSeqAIJ(PETSC_COMM_SELF,ySize*nz,ny*nz,9,PETSC_NULL,&DAdyseq);CHKERRQ(ierr);
    ierr  = MatCreateSeqAIJ(PETSC_COMM_SELF,ySize*nz,ny*nz,9,PETSC_NULL,&DAdzseq);CHKERRQ(ierr);
  
        ind = 0;
        for(i=0;i<ySize;i++){
              for (j=0;j<nz;j++){
             
                  p = 2*i+1;
                  q = 2*j+1;
   
                  vr = *(velx0+(i+yStart)*nz+j); 

                  for(k=0;k<9;k++){
                     p   = 2*i+1+*(yshift+k);
                     q   = 2*j+1+*(zshift+k);
                     *(ye+k) = *(y0+p*(2*nz+1)+q);
                     *(ze+k) = *(z0+p*(2*nz+1)+q); 
                     *(iy+k) = (PetscInt)floor(*(ye+k)/gsizey);
                     *(iz+k) = (PetscInt)floor(*(ze+k)/gsizez);   

                   
                     if(cdinfo->period->y==0){
                      if(*(iy+k)>=ny){*(iy+k)=ny-1;}
                      if(*(iy+k)<0){*(iy+k)=0;}
                     }else{
                      if(*(iy+k)>=ny){*(iy+k)=*(iy+k)-ny;}
                      if(*(iy+k)<0){*(iy+k)=*(iy+k)+ny;}
                     }

                     if(cdinfo->period->z==0){
                         if(*(iz+k)>=nz){*(iz+k)=nz-1;}                      
                         if(*(iz+k)<0){*(iz+k)=0;}
                     }else{
                         if(*(iz+k)>=nz){*(iz+k)=*(iz+k)-nz;}                      
                         if(*(iz+k)<0){*(iz+k)=*(iz+k)+nz;}
                      }
     

                     *(col+k) = *(iy+k)*nz+*(iz+k);
                     for(l=0;l<k;l++){ if(*(col+l)==*(col+k)) {*(col+k) =-1;}}
                  }
  
                  valsum = 0; 
                  for(k=0;k<9;k++){
                      if(*(col+k)>0){
                         py = *(ye+4)-(*(iy+k)*gsizey+hgsizey);
                         pz = *(ze+4)-(*(iz+k)*gsizez+hgsizez);
        
                         *(val+k)=fabs(gsizey-fabs(py))*fabs(gsizez-fabs(pz));
                         if(py>0){*(dAdyval+k) =  fabs(gsizez-fabs(pz));}
                             else{*(dAdyval+k) = -fabs(gsizez-fabs(pz));} 
                         if(pz>0){*(dAdzval+k) =  fabs(gsizey-fabs(py));}
                             else{*(dAdzval+k) = -fabs(gsizey-fabs(py));} 
                         valsum = valsum+*(val+k)*(*(velx0+*(iy+k)*nz+*(iz+k)));
                      }else{
                         *(val+k) = 0;
                      }
                  }
                  for(k=0;k<9;k++){*(val+k)= *(val+k)/valsum*vr; if(isnan(*(val+k))){*(val+k)=1;} }
            
                  ierr = MatSetValues(Aseq,1,&ind,9,col,val,ADD_VALUES);CHKERRQ(ierr);
                  ierr = MatSetValues(DAdyseq,1,&ind,9,col,dAdyval,ADD_VALUES);CHKERRQ(ierr);
                  ierr = MatSetValues(DAdzseq,1,&ind,9,col,dAdzval,ADD_VALUES);CHKERRQ(ierr);
                  ind = ind+1;
             }
        }

        ierr = MatAssemblyBegin(Aseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(Aseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
        ierr = MatAssemblyBegin(DAdyseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(DAdyseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyBegin(DAdzseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(DAdzseq,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

        ierr =  MatMerge(MPI_COMM_WORLD,Aseq,ySize*nz,MAT_INITIAL_MATRIX,&A);CHKERRQ(ierr);
        ierr =  MatMerge(MPI_COMM_WORLD,DAdyseq,ySize*nz,MAT_INITIAL_MATRIX,&DAdy);CHKERRQ(ierr);
        ierr =  MatMerge(MPI_COMM_WORLD,DAdzseq,ySize*nz,MAT_INITIAL_MATRIX,&DAdz);CHKERRQ(ierr);
       
        *MapA = A; 
        *dAdy = DAdy;
        *dAdz = DAdz;
   

////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////


        ierr = VecDestroy(x0vec);CHKERRQ(ierr);
        ierr = VecDestroy(y0vec);CHKERRQ(ierr);
        ierr = VecDestroy(z0vec);CHKERRQ(ierr);
        ierr = VecDestroy(vx0vec);CHKERRQ(ierr);
        free(ye);free(ze);free(val); 
        free(yshift); free(zshift);
        free(iy);free(iz);free(col);
        
  
return 0;

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
int dMarkovMap2d(Vec sol, PetscInt ny, PetscInt nz, CoordInfo *cdinfo,Mat *dydv, Mat *dzdv)
{

   Vec            temp, x0vec, y0vec, z0vec;
   Mat            Dydv, Dzdv;
   PetscErrorCode ierr;
   PetscScalar    gsizey, gsizez;
   PetscInt       i, j, n, np, ind, yStart, yEnd, ySize;
   PetscScalar    *x0, *y0, *z0;

    gsizey = (cdinfo->r->y-cdinfo->l->y)/(1.0*ny);
    gsizez = (cdinfo->r->z-cdinfo->l->z)/(1.0*nz);
   n      = 500;

   ierr = VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,ny,&temp);
   ierr = VecGetOwnershipRange(temp,&yStart,&yEnd);    
   ySize = yEnd-yStart;


    np      = ySize*nz;
    x0      = malloc(np*sizeof(PetscScalar));
    y0      = malloc(np*sizeof(PetscScalar));
    z0      = malloc(np*sizeof(PetscScalar));

    ind = 0;
    for(i=yStart;i<yEnd;i++){
         for (j=0;j< nz;j++){
                *(x0+ind) = 0;
                *(y0+ind) = (i+0.5)*gsizey;
                *(z0+ind) = (j+0.5)*gsizez;      
                ind = ind+1;
         }
    }
    // use seq vector
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)x0,&x0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)y0,&y0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)z0,&z0vec);CHKERRQ(ierr);


   StreamlineMapr(&x0vec, &y0vec,   &z0vec,   sol, cdinfo, &Dydv, &Dzdv,32,  n ,10);
  
   
   *dydv = Dydv;
   *dzdv = Dzdv;  

return 0;

}

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//   This function returns  df/dv which f is now designed to be a twist map
// 
//   Tzu-Chen Liang 9-17-2007

int dMarkovMap2df(Vec sol, PetscInt ny, PetscInt nz, CoordInfo *cdinfo,Vec *dfdv)
{

   Vec            temp, x0vec, y0vec, z0vec;
   Vec            py0vec,pz0vec;
   Mat            Dydv, Dzdv;
   Vec            Dfdv;
   PetscErrorCode ierr;
   PetscScalar    gsizey, gsizez;
   PetscInt       i, j, n, np, ind, yStart, yEnd, ySize;
   PetscScalar    *x0, *y0, *z0, yind,zind;
   PetscScalar    *xe, *ye, *ze;
   Vec            xevec, yevec, zevec; 
   PetscScalar    yc,zc,alpha,r,th;
   PetscInt       yg,zg;  
   PetscInt       solmaxp;
   PetscScalar    solmax;

   ierr = VecMax(sol,&solmaxp,&solmax);CHKERRQ(ierr);
   
   gsizey = (cdinfo->r->y-cdinfo->l->y)/(1.0*ny);
   gsizez = (cdinfo->r->z-cdinfo->l->z)/(1.0*nz);
   n      = 2000;

   ierr = VecCreateMPI(MPI_COMM_WORLD,PETSC_DECIDE,ny,&temp);
   ierr = VecGetOwnershipRange(temp,&yStart,&yEnd);    
   ySize = yEnd-yStart;
   //ierr = VecDestroy(temp);

    np      = ySize*nz;
    x0      = malloc(np*sizeof(PetscScalar));
    y0      = malloc(np*sizeof(PetscScalar));
    z0      = malloc(np*sizeof(PetscScalar));

    ind = 0;
    for(i=yStart;i<yEnd;i++){
         for (j=0;j< nz;j++){
                *(x0+ind) = 0;
                *(y0+ind) = (i+0.5)*gsizey;
                *(z0+ind) = (j+0.5)*gsizez;      
                ind = ind+1;
         }
    }
    // use seq vector
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)x0,&x0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)y0,&y0vec);CHKERRQ(ierr);
    ierr  = VecCreateSeqWithArray(PETSC_COMM_SELF,np,(PetscScalar*)z0,&z0vec);CHKERRQ(ierr);

   
  
   ierr = VecDuplicate(x0vec,&xevec);CHKERRQ(ierr); 
   ierr = VecDuplicate(y0vec,&yevec);CHKERRQ(ierr);
   ierr = VecDuplicate(z0vec,&zevec);CHKERRQ(ierr);
   ierr = VecCopy(x0vec,xevec);CHKERRQ(ierr);
   ierr = VecCopy(y0vec,yevec);CHKERRQ(ierr);
   ierr = VecCopy(z0vec,zevec);CHKERRQ(ierr);
 
   // dt =0.01/solmax
   StreamlineMapr(&xevec, &yevec,   &zevec,   sol, cdinfo, &Dydv, &Dzdv,4*0.01/solmax,  n ,10);
  // StreamlineMapr(&x0vec, &y0vec,   &z0vec,   sol, cdinfo, &Dydv, &Dzdv,32,  n ,10);
  // SparseSave(Dydv,"Dydv");
  // SparseSave(Dzdv,"Dzdv");
 
   ierr = VecGetArray(xevec,&xe);CHKERRQ(ierr);
   ierr = VecGetArray(yevec,&ye);CHKERRQ(ierr);
   ierr = VecGetArray(zevec,&ze);CHKERRQ(ierr);

   yc = 0.35;
   zc = 0.35;
   alpha = M_PI/4;
   ind=0;
   // a ccw gradient
   for(i=yStart;i<yEnd;i++){
         for (j=0;j< nz;j++){

// standard map

        //     yind = *(ye+ind)- (*(y0+ind)+*(z0+ind)+0.3*sin(2*M_PI*(*(y0+ind)))); 
        //     zind = *(ze+ind)- (*(z0+ind)+0.3*sin(2*M_PI*(*(y0+ind))));
          




//  finite angle rotation

              r = sqrt(pow(*(y0+ind)-yc, 2)  + pow(*(z0+ind)-zc,2));
              th = atan2( *(z0+ind) - zc,*(y0+ind) - yc)+alpha;                 
              //if(th>2*M_PI) {th =th-2*M_PI; }
              if (r<0.35){
                yind = (*(ye+ind)-yc-r*cos(th));
                zind = (*(ze+ind)-zc-r*sin(th));
		  }else{
                yind = *(ye+ind)-*(y0+ind);
                zind = *(ze+ind)-*(z0+ind); 
                }



         //y[0 1], z[0 1]

/*
        yg = floor(*(y0+ind)*2);
        zg = floor(*(z0+ind)*2);
        if(yg>1){yg = 1;}
        if(zg>1){zg = 1;}


        if((yg==0)&&(zg==0)){//right 
           yind = *(ye+ind)-*(y0+ind)+0.5;
           zind = *(ze+ind)-*(z0+ind);
        }

        if((yg==1)&&(zg==0)) {//up
           yind = *(ye+ind)-*(y0+ind);
           zind = *(ze+ind)-*(z0+ind)+0.5;
        }
        if((yg==1)&&(zg==1)){//left
           yind = *(ye+ind)-*(y0+ind)-0.5;
           zind = *(ze+ind)-*(z0+ind);
        }

        if((yg==0)&&(zg==1)){//down
           yind = *(ye+ind)-*(y0+ind);
           zind = *(ze+ind)-*(z0+ind)-0.5;
        } 
*/
  /*    
        yg = floor(*(y0+ind)*4);
        zg = floor(*(z0+ind)*4);
        if(yg>3){yg = 3;}
        if(zg>3){zg = 3;}


        if((zg==0)&&(yg==(1||2||3))){//left 
           yind = *(ye+ind)-*(y0+ind)-0.25;
           zind = *(ze+ind)-*(z0+ind);
        }

        if((yg==0)&&(zg==(0||1||2))) {//up
           yind = *(ye+ind)-*(y0+ind);
           zind = *(ze+ind)-*(z0+ind)+0.25;
        }
        if((yg==2)&&(zg==(1||2))){//up
           yind =15*( *(ye+ind)-*(y0+ind));
           zind =15*( *(ze+ind)-*(z0+ind))+0.25;
        }

        if((yg==1)&&(zg==(2||3))){//down
           yind =15*(*(ye+ind)-*(y0+ind));
           zind =15*(*(ze+ind)-*(z0+ind)-0.25);
        }
        if((yg==3)&&(zg==(1||2||3))){//down
           yind = *(ye+ind)-*(y0+ind);
           zind = *(ze+ind)-*(z0+ind)-0.25;
        }

        if(((yg==0)&&(zg==3))||((yg==1)&&(zg==1))||((yg==2)&&(zg=3))){//right
           yind = *(ye+ind)-*(y0+ind)+0.25;
           zind = *(ze+ind)-*(z0+ind);
        } 
       

   */  








    
    /*
       if(*(z0+ind)<1.0/3) {yind = -1; zind = 0;}
           elseif(*(z0+ind)>2.0/3) {yind =  1; zind = 0;}
           elseif(*(y0+ind)<1.0/3) {yind =  0; zind = 1;}
     elseif(*(y0+ind)<1.0&&*(y0+ind)>2.0/3) {yind =  0; zind = 1;}
     elseif(*(y0+ind)<5.0/3&&*(y0+ind)>4.0/3) {yind =  0; zind = 1;}        
     else{yind =  0; zind = -1;}      
    */

          /*   
             if(*(y0+ind)<1 && *(z0+ind)<0.5){
                     yind = *(ye+ind)-*(y0+ind)-0.5;
                     zind = *(ze+ind)-*(z0+ind)-0.5;
                     }
             if(*(y0+ind)<1 && *(z0+ind)>0.5){
                     yind = *(ye+ind)-*(y0+ind)-0.5;
                     zind = *(ze+ind)-*(z0+ind)+0.5;
                     }
             if(*(y0+ind)>1 && *(z0+ind)<0.5){
                     yind = *(ye+ind)-*(y0+ind)+0.5;
                     zind = *(ze+ind)-*(z0+ind)-0.5;
                     }
             if(*(y0+ind)>1 && *(z0+ind)>0.5){
                     yind = *(ye+ind)-*(y0+ind)+0.5;
                     zind = *(ze+ind)-*(z0+ind)+0.5;
                     }

           */ 






         /*             
              r = sqrt(pow(*(y0+ind)-yc, 2)  + pow(*(z0+ind)-zc,2));
              th = atan2( *(z0+ind) - zc,*(y0+ind) - yc)+0.3;                 
              //if(th>2*M_PI) {th =th-2*M_PI; }
              if (r<2){
                yind = (*(ye+ind)-yc-r*cos(th));
                zind = (*(ze+ind)-zc-r*sin(th));
		  }else{
                yind = 0;
                zind = 0; 
                }
         */

		 // fprintf(stdout,"yindzind  = %f %f \n", yind , zind );
               //   yind = (*(z0+ind)-0.5);
               //   zind =-(*(y0+ind)-0.5);
                *(y0+ind)=yind;  
                *(z0+ind)=zind;                    
                ind = ind+1;
         }
    }

  // use y0,z0 as output
  ierr  = VecCreateMPIWithArray(PETSC_COMM_WORLD,np,PETSC_DECIDE,y0,&py0vec);CHKERRQ(ierr);
  ierr  = VecCreateMPIWithArray(PETSC_COMM_WORLD,np,PETSC_DECIDE,z0,&pz0vec);CHKERRQ(ierr);



  ierr  =  MatGetVecs(Dydv,&Dfdv,PETSC_NULL);CHKERRQ(ierr);
  ierr  =  MatMultTranspose(Dydv,py0vec,Dfdv);CHKERRQ(ierr);
  ierr  =  MatMultTransposeAdd(Dzdv,pz0vec,Dfdv,Dfdv);CHKERRQ(ierr);
 


   *dfdv = Dfdv;

ierr =VecDestroy(x0vec);    
ierr =VecDestroy(y0vec);
ierr =VecDestroy(z0vec);
free(x0);free(y0);free(z0);

ierr =VecDestroy(xevec);
ierr =VecDestroy(yevec);
ierr =VecDestroy(zevec);
ierr =MatDestroy(Dydv);
ierr =MatDestroy(Dzdv); 

return 0;

}


