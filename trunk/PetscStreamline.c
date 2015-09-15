#include "petscksp.h"
#include "petscda.h"
#include "slepceps.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "MeshStruc.h"
#include "mpi.h"
#include "PetscReceive.h"
 
#define DIMx 0
#define DIMy 1
#define DIMz 2
#define DIRx 0
#define DIRy 1
#define DIRz 2

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

int loc2ijk(PetscInt vdir, PetscInt dir, PetscInt period, PetscScalar gsize, PetscScalar x, PetscInt *indi, PetscScalar *rem, PetscInt gridspam)
{
  PetscScalar temp,remi;
  PetscInt    ijk;
  

  if(vdir==dir){
           
      temp  =  x + gsize;
      ijk   =  floor(temp/gsize);
      remi  =  temp - ijk*gsize;
      if (period == 0) { ijk = ijk -1;}  
      if (period == 1 && ijk ==gridspam) {ijk =0 ; }   
 
  }else{
      temp  = x + 0.5*gsize;
      ijk   = floor(temp/gsize);
      remi  = temp -ijk*gsize;
  }
      ijk = ijk -1; 
      
      remi = remi / gsize; // We normalize remainder here!
   
       

      *indi = ijk;
      *rem  = remi;
return(0);
}

///////////////////////////////////////////////////////////////////////////////////////
int gloc2ijk(PetscInt3 *period, PetscScalar3 *gsize,PetscScalar x,PetscScalar y,PetscScalar z, PetscInt *ind, PetscScalar *rem, PetscInt33* gridspam)
{

   loc2ijk(DIMx,DIRx,period->x,gsize->x,x,ind  ,rem  ,gridspam->u.x);
   loc2ijk(DIMx,DIRy,period->x,gsize->x,y,ind+1,rem+1,gridspam->u.y);
   loc2ijk(DIMx,DIRz,period->x,gsize->x,z,ind+2,rem+2,gridspam->u.z);
   loc2ijk(DIMy,DIRx,period->y,gsize->y,x,ind+3,rem+3,gridspam->v.x);
   loc2ijk(DIMy,DIRy,period->y,gsize->y,y,ind+4,rem+4,gridspam->v.y);
   loc2ijk(DIMy,DIRz,period->y,gsize->y,z,ind+5,rem+5,gridspam->v.z);
   loc2ijk(DIMz,DIRx,period->z,gsize->z,x,ind+6,rem+6,gridspam->w.x);
   loc2ijk(DIMz,DIRy,period->z,gsize->z,y,ind+7,rem+7,gridspam->w.y);
   loc2ijk(DIMz,DIRz,period->z,gsize->z,z,ind+8,rem+8,gridspam->w.z);



  return(0);
}

///////////////////////////////////////////////////////////////////////////////////////
int ijk2ind(PetscInt *ijk,PetscInt3 *gridspam, PetscInt3 *period ,PetscInt *ind)
{
  PetscInt i=0,j=0; 
  

     for(j=0;j<24;j=j+3){

           if(period->x==1) {
             if(ijk[j]<0)             {ijk[j]   = ijk[j]+gridspam->x; } 
             if(ijk[j]>=gridspam->x)   {ijk[j]   = ijk[j]-gridspam->x; }
           }
           if(period->y==1) {
             if(ijk[j+1]<0)           {ijk[j+1] = ijk[j+1]+gridspam->y; } 
             if(ijk[j+1]>=gridspam->y) {ijk[j+1] = ijk[j+1]-gridspam->y; }
           }
           if(period->z==1) {
             if(ijk[j+2]<0)           {ijk[j+2] = ijk[j+2]+gridspam->z; } 
             if(ijk[j+2]>=gridspam->z) {ijk[j+2] = ijk[j+2]-gridspam->z; }
           }

          
           ind[i] = ijk[j]%(gridspam->x) +0; // x dir is periodical...
           ind[i] = ind[i] + ijk[j+1]*gridspam->x;
           ind[i] = ind[i] + ijk[j+2]*gridspam->x*gridspam->y;  
           
            
         if(period->y==0){
                if((ijk[j+1]<0)||(ijk[j+1]>=gridspam->y)) {ind[i]=-1;}
           }
         if(period->z==0){
                if((ijk[j+2]<0)||(ijk[j+2]>=gridspam->z)) {ind[i]=-1;}
           }


          //if(*period==0){
           if(ijk[j]<0||ijk[j+1]<0||ijk[j+2]<0){ind[i]=-1;}
           if(ijk[j+1]>=gridspam->y||ijk[j+2]>=gridspam->z){ind[i]=-1;}
         // }
           
         

           i      =  i + 1; 
      
            

  }
 
  return(0);
}
///////////////////////////////////////////////////////////////////////////////////////
int gijk2ind(PetscInt *ijk,PetscInt33 *gridspam,PetscInt3 *period, PetscInt *ind)
{

   ijk2ind(ijk   ,&(gridspam->u),period  ,ind);
   ijk2ind(ijk+24,&(gridspam->v),period  ,ind+8);
   ijk2ind(ijk+48,&(gridspam->w),period  ,ind+16);
   
 return(0);
}
///////////////////////////////////////////////////////////////////////////////////////
int cornerfind(PetscInt *ijk, PetscInt *cornerijk)
{
  PetscInt j;
  PetscInt iarray[24]={0,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,1,0,1,1,1};
  
     PetscMemcpy(cornerijk   ,ijk      , 3*sizeof(PetscInt));
     PetscMemcpy(cornerijk+3 ,cornerijk, 3*sizeof(PetscInt));
     PetscMemcpy(cornerijk+6 ,cornerijk, 6*sizeof(PetscInt));
     PetscMemcpy(cornerijk+12,cornerijk,12*sizeof(PetscInt));
     for(j=0;j<24;j++){cornerijk[j] = cornerijk[j] + iarray[j];}

  return(0);

}

///////////////////////////////////////////////////////////////////////////////////////
int gcornerfind(PetscInt *ijk, PetscInt *cornerijk)
{
   cornerfind(ijk   ,cornerijk);
   cornerfind(ijk+3 ,cornerijk+24);
   cornerfind(ijk+6 ,cornerijk+48);

   return(0);
}
///////////////////////////////////////////////////////////////////////////////////////
int lind2gind(PetscInt *ind,PetscInt dir,PetscInt gridshift,PetscInt zeroind)
{
   PetscInt i;
   
   for(i=0;i<8;i++){
    
    if (ind[i]>=0) {ind[i] = ind[i] + gridshift;}
    else{ind[i]= zeroind+dir*8+i;} 
   
   }
    

   return(0);
}
///////////////////////////////////////////////////////////////////////////////////////
int glind2gind(PetscInt *ind,PetscInt4 *nofgrid, PetscInt zeroind)
{
   lind2gind(ind   ,DIMx,0,zeroind); 
   lind2gind(ind+8 ,DIMy,nofgrid->x,zeroind); 
   lind2gind(ind+16,DIMz,nofgrid->x+nofgrid->y,zeroind); 
   return(0);
}
///////////////////////////////////////////////////////////////////////////////////////
int lininterp3d(PetscScalar *uval,PetscScalar *x,PetscScalar *u,PetscScalar *w)
{
   PetscInt     i;
   PetscScalar  xn[3];
  
    

   for(i=0;i<3;i++){xn[i] = 1- x[i];}
   w[0] = xn[0]*xn[1]*xn[2]; 
   w[1] = xn[0]*xn[1]* x[2];
   w[2] = xn[0]* x[1]*xn[2];
   w[3] = xn[0]* x[1]* x[2];
   w[4] =  x[0]*xn[1]*xn[2];
   w[5] =  x[0]*xn[1]* x[2];
   w[6] =  x[0]* x[1]*xn[2];
   w[7] =  x[0]* x[1]* x[2];

   *u = 0;
   for(i=0;i<8;i++) {*u = *u + w[i]*uval[i]; }
 
    
 
   return(0);
}
///////////////////////////////////////////////////////////////////////////////////////
int glininterp3d(PetscScalar *uval,PetscScalar *rem,PetscScalar3 *u,PetscScalar *w)
{
  

   lininterp3d(uval   ,rem  ,&(u->x)  ,w);
   lininterp3d(uval+8 ,rem+3,&(u->y)  ,w+8);
   lininterp3d(uval+16,rem+6,&(u->z)  ,w+16);
  

   return(0);
}

///////////////////////////////////////////////////////////////////////////////////////
int veval(PetscScalar x, PetscScalar y, PetscScalar z, Vec *mysol, CoordInfo *cdinfo, PetscScalar3 *vel, PetscInt *ind, PetscScalar *weight)
{

   PetscInt      ijk[9],cornerijk[72];
   PetscScalar   uval[9];
   PetscScalar   remainder[9];
   PetscScalar3  uvw;
   PetscInt      zeroind;   
   
   zeroind = cdinfo->nofgrid->x+cdinfo->nofgrid->y+cdinfo->nofgrid->z;
   gloc2ijk(cdinfo->period,cdinfo->gsize,x,y,z,ijk,remainder,cdinfo->gridspam);
   gcornerfind(ijk,cornerijk);
   gijk2ind(cornerijk,cdinfo->gridspam,cdinfo->period,ind);
   glind2gind(ind,cdinfo->nofgrid,zeroind);  
   VecGetValues(*mysol,24,ind,uval);
   glininterp3d(uval,remainder,&uvw,weight);

   
   
 

    *vel = uvw;
 


  return(0);
}
///////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int RK2(PetscScalar *x, PetscScalar *y, PetscScalar *z, Vec *sol, CoordInfo *cdinfo, Mat *rowmat, PetscInt j,PetscScalar dt)
{
     PetscErrorCode ierr;
     PetscInt       k;
     PetscInt       ind[24];
     PetscScalar    weight[24];
     PetscScalar3   k1;
     PetscScalar3   vel;
//printf(stdout,"x=%f, y=%f, z=%f \n",*x,*y,*z);
          veval(*x, *y, *z, sol, cdinfo, &vel, ind, weight);             
          for(k=0;k<24;k++){weight[k] = dt*weight[k];}        

         if(rowmat){ 
          ierr = MatSetValues(*rowmat,1,&j,24,ind,weight,ADD_VALUES);CHKERRQ(ierr);
          }

          k1.x =  dt * vel.x;
          k1.y =  dt * vel.y;
          k1.z =  dt * vel.z;

          *x = *x + k1.x ;
          *y = *y + k1.y ;
          *z = *z + k1.z ;
        
         
          if(*y>cdinfo->r->y){
            if(cdinfo->period->y==0){*y=cdinfo->r->y;}
                                else{*y=*y-cdinfo->r->y;} 
          }

          if(*y<cdinfo->l->y){
            if(cdinfo->period->y==0){*y=cdinfo->l->y;}
                                else{*y=*y+cdinfo->r->y;}
          }
        
          if(*z>cdinfo->r->z){
            if(cdinfo->period->z==0){*z=cdinfo->r->z;}
                                else{*z=*z-cdinfo->r->z;}
          }

          if(*z<cdinfo->l->z){
              if(cdinfo->period->z==0){*z=cdinfo->l->z;}
                                else{*z=*z+cdinfo->r->z;} 
          }
         

return 0;

}


/////////////////////////////////////////////////////////////////////////////////
int RK4(PetscScalar *x, PetscScalar *y, PetscScalar *z, Vec *sol, CoordInfo *cdinfo, Mat *rowmat, PetscInt j,PetscScalar dt)
{
     PetscErrorCode ierr;
     PetscInt       k;
     PetscInt       ind[24];
     PetscScalar    weight[24];
     PetscScalar3   k1,k2,k3,k4;
     PetscScalar3   vel;

        veval(*x, *y, *z, sol, cdinfo, &vel, ind, weight);             
          for(k=0;k<24;k++){weight[k] = dt/6*weight[k];}                  
          ierr = MatSetValues(*rowmat,1,&j,24,ind,weight,ADD_VALUES);CHKERRQ(ierr);        
          k1.x =  dt * vel.x;
          k1.y =  dt * vel.y;
          k1.z =  dt * vel.z;

        veval(*x+0.5*k1.x, *y+0.5*k1.y, *z+0.5*k1.z, sol, cdinfo, &vel, ind, weight);
          for(k=0;k<24;k++){weight[k] = dt/3*weight[k];}  
          ierr = MatSetValues(*rowmat,1,&j,24,ind,weight,ADD_VALUES);CHKERRQ(ierr);
          k2.x =  dt * vel.x;
          k2.y =  dt * vel.y;
          k2.z =  dt * vel.z;

        veval(*x+0.5*k2.x, *y+0.5*k2.y, *z+0.5*k2.z, sol, cdinfo, &vel, ind, weight);
          for(k=0;k<24;k++){weight[k] = dt/3*weight[k];}  
          ierr = MatSetValues(*rowmat,1,&j,24,ind,weight,ADD_VALUES);CHKERRQ(ierr);
          k3.x =  dt * vel.x;
          k3.y =  dt * vel.y;
          k3.z =  dt * vel.z;

        veval(*x + k3.x, *y + k3.y, *z + k3.z, sol, cdinfo, &vel, ind, weight);
          for(k=0;k<24;k++){weight[k] = dt/6*weight[k];}  
          ierr = MatSetValues(*rowmat,1,&j,24,ind,weight,ADD_VALUES);CHKERRQ(ierr);
          k4.x =  dt * vel.x;
          k4.y =  dt * vel.y;
          k4.z =  dt * vel.z;

 
          *x = *x + (k1.x)/6 + (k2.x)/3 + (k3.x)/3 + (k4.x)/6;
          *y = *y + (k1.y)/6 + (k2.y)/3 + (k3.y)/3 + (k4.y)/6;
          *z = *z + (k1.z)/6 + (k2.z)/3 + (k3.z)/3 + (k4.z)/6;
    

 return 0;

}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int StreamlineMap(Vec x0, Vec y0, Vec z0, Vec v, CoordInfo *cdinfo, Mat *Dydv, Mat *Dzdv,PetscScalar dt, PetscInt n ,PetscInt trycountmax, PetscScalar **xstream, PetscScalar **ystream, PetscScalar **zstream) 
{
   PetscErrorCode ierr;
   PetscScalar    endx; 
   PetscInt       solsize ,trycount = 0;
   PetscInt       i,j,k,Istart,Iend,Lsize,Ivstart,Ivend;
   Vec            dydx, dzdx, Ivec1, Ivec2, Ivec3;
   Mat            rowmat,weightmat, Dx, Dy, Dz;
   PetscInt     ind[24];
   PetscScalar  weight[24];
   PetscScalar3 vel;
 
   ierr = VecGetOwnershipRange(x0,&Istart,&Iend);CHKERRQ(ierr); 
   Lsize = Iend - Istart; 


   PetscScalar xlast,ylast,zlast; 
   PetscScalar *xarray, *yarray,*zarray;  
   PetscScalar steplen;

    PetscScalar *x,*y,*z;
    x = malloc(Lsize*n*sizeof(PetscScalar));
    y = malloc(Lsize*n*sizeof(PetscScalar));
    z = malloc(Lsize*n*sizeof(PetscScalar));

   



   solsize   = cdinfo->nofgrid->x + cdinfo->nofgrid->y + cdinfo->nofgrid->z + 24;
   endx      = cdinfo->r->x;
   // add this to prevent the last grid problem..
   endx      = endx - cdinfo->gsize->x;
  
   ierr = VecGetArray(x0,&xarray);CHKERRQ(ierr);
   ierr = VecGetArray(y0,&yarray);CHKERRQ(ierr);
   ierr = VecGetArray(z0,&zarray);CHKERRQ(ierr);
  

   ierr = VecCreateMPI(PETSC_COMM_WORLD,Lsize,PETSC_DETERMINE,&dydx);CHKERRQ(ierr); 
   ierr = VecCreateMPI(PETSC_COMM_WORLD,Lsize,PETSC_DETERMINE,&dzdx);CHKERRQ(ierr); 

   ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,solsize,&Ivec1);CHKERRQ(ierr); 
   ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,solsize,&Ivec2);CHKERRQ(ierr);
   ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,solsize,&Ivec3);CHKERRQ(ierr);

   ierr = VecGetOwnershipRange(Ivec1,&Ivstart,&Ivend);CHKERRQ(ierr);
   for(i=0;i<cdinfo->nofgrid->x;i++){
         if(i>=Ivstart&&i<Ivend){ierr = VecSetValue(Ivec1,i,1, INSERT_VALUES);}
                   } 
   for(i=cdinfo->nofgrid->x;i<cdinfo->nofgrid->x+cdinfo->nofgrid->y;i++){
         if(i>=Ivstart&&i<Ivend){ierr = VecSetValue(Ivec2,i,1, INSERT_VALUES);}
                   } 
   for(i=cdinfo->nofgrid->x+cdinfo->nofgrid->y;i<solsize-24;i++) {
         if(i>=Ivstart&&i<Ivend){ierr = VecSetValue(Ivec3,i,1, INSERT_VALUES);}
                   } 
  

   ierr = VecAssemblyBegin(Ivec1);CHKERRQ(ierr); 
   ierr = VecAssemblyEnd(Ivec1);CHKERRQ(ierr); 
   ierr = VecAssemblyBegin(Ivec2);CHKERRQ(ierr); 
   ierr = VecAssemblyEnd(Ivec2);CHKERRQ(ierr); 
   ierr = VecAssemblyBegin(Ivec3);CHKERRQ(ierr); 
   ierr = VecAssemblyEnd(Ivec3);CHKERRQ(ierr); 


// Now let's begin the streamline    

   ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,Lsize,solsize,(PetscInt)(40*cdinfo->gridspam->u.x),PETSC_NULL,&rowmat);CHKERRQ(ierr);

   for(j=0;j<Lsize;j++){
       *(x+j*n) = *(xarray+j);
       *(y+j*n) = *(yarray+j);
       *(z+j*n) = *(zarray+j);
       xlast   = *(x+j*n);
       ylast   = *(y+j*n);
       zlast   = *(z+j*n);
       i       = 0;

       // test boundary
       //if(pow(*(y+j*n)-cdinfo->l->y,2)<1e-4) {*(x+j*n+i)=endx;}
       //if(pow(*(y+j*n)-cdinfo->r->y,2)<1e-4) {*(x+j*n+i)=endx;}
       //if(pow(*(z+j*n)-cdinfo->l->z,2)<1e-4) {*(x+j*n+i)=endx;}
       //if(pow(*(z+j*n)-cdinfo->l->z,2)<1e-4) {*(x+j*n+i)=endx;}    
  
       while(*(x+j*n+i)<endx && i<n-1){      
             RK4(&xlast, &ylast, &zlast, &v, cdinfo, &rowmat, j,dt);    
             steplen = sqrt(pow(*(x+j*n+i)-xlast,2) + pow(*(y+j*n+i)-ylast,2)+pow(*(z+j*n+i)-zlast,2));

             if(steplen >=2*(endx*1.0)/n||trycount>trycountmax){
                    *(x+j*n+i+1) = xlast;
                    *(y+j*n+i+1) = ylast;
                    *(z+j*n+i+1) = zlast;
                    trycount = 0;
                    i++;
             }else{
                    trycount++; 
             }

       }
 
      for(k=i;k<n-1;k++){
  
            *(x+j*n+k+1)=*(x+j*n+k);
            *(y+j*n+k+1)=*(y+j*n+k);
            *(z+j*n+k+1)=*(z+j*n+k);

         }
     
     veval(*(x+j*n+n-1),*(y+j*n+n-1), *(z+j*n+n-1), &v,cdinfo, &vel, ind, weight);
     
    
     
  
     VecSetValue(dydx,Istart+j,(PetscScalar)vel.y/vel.x, INSERT_VALUES); 
     VecSetValue(dzdx,Istart+j,(PetscScalar)vel.z/vel.x, INSERT_VALUES); 
     // Replace the locations by the end points 
     VecSetValue(x0,Istart+j,*(x+j*n+n-1), INSERT_VALUES); 
     VecSetValue(y0,Istart+j,*(y+j*n+n-1), INSERT_VALUES); 
     VecSetValue(z0,Istart+j,*(z+j*n+n-1), INSERT_VALUES); 



} 




      ierr = VecRestoreArray(x0,&xarray);
      ierr = VecRestoreArray(y0,&yarray); 
      ierr = VecRestoreArray(z0,&zarray);


      ierr = VecAssemblyBegin(dydx);
      ierr = VecAssemblyEnd(dydx);
      ierr = VecAssemblyBegin(dzdx);
      ierr = VecAssemblyEnd(dzdx);
      ierr = VecAssemblyBegin(x0);
      ierr = VecAssemblyEnd(x0);
      ierr = VecAssemblyBegin(y0);
      ierr = VecAssemblyEnd(y0);
      ierr = VecAssemblyBegin(z0);
      ierr = VecAssemblyEnd(z0);

      ierr = MatAssemblyBegin(rowmat,MAT_FINAL_ASSEMBLY);
      ierr = MatAssemblyEnd(rowmat,MAT_FINAL_ASSEMBLY);

  
      MatMerge(PETSC_COMM_WORLD,rowmat,PETSC_DECIDE,MAT_INITIAL_MATRIX,&weightmat);

      MatDuplicate(weightmat,MAT_COPY_VALUES,&Dx); 
      MatDuplicate(weightmat,MAT_COPY_VALUES,&Dy);
      MatDuplicate(weightmat,MAT_COPY_VALUES,&Dz);
   
  
      // Dy = Dy + dy/dx * Dx
      MatDiagonalScale(Dx,dydx,Ivec1);
      MatDiagonalScale(Dy,PETSC_NULL,Ivec2); 
      MatAXPY(Dy,1,Dx,SAME_NONZERO_PATTERN);
      // Dz = Dz + dz/dx * Dx
      MatDuplicate(weightmat,MAT_COPY_VALUES,&Dx); 
      MatDiagonalScale(Dx,dzdx,Ivec1);
      MatDiagonalScale(Dz,PETSC_NULL,Ivec3);  
      MatAXPY(Dz,1,Dx, SAME_NONZERO_PATTERN);


      *xstream = (PetscScalar*)x;
      *ystream = (PetscScalar*)y;
      *zstream = (PetscScalar*)z;
      




      *Dydv = Dy;
      *Dzdv = Dz;

   ierr = VecDestroy(dydx);CHKERRQ(ierr);
   ierr = VecDestroy(dzdx);CHKERRQ(ierr);

   ierr = VecDestroy(Ivec1);CHKERRQ(ierr);
   ierr = VecDestroy(Ivec2);CHKERRQ(ierr);
   ierr = VecDestroy(Ivec3);CHKERRQ(ierr);
   //ierr = MatDestroy(weightmat);CHKERRQ(ierr);




return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int StreamlineMapr(Vec *x0, Vec *y0, Vec *z0, Vec v, CoordInfo *cdinfo, Mat *Dydv, Mat *Dzdv,PetscScalar dt, PetscInt n ,PetscInt trycountmax) 
{
   PetscErrorCode ierr;
   PetscScalar    endx,xc,yc,zc; 
   PetscInt       solsize ,trycount = 0;
   PetscInt       i,j,Istart,Iend,Lsize,Ivstart,Ivend;
   Vec            dydx, dzdx, Ivec1, Ivec2, Ivec3;
   Mat            rowmat,weightmat, Dx, Dy, Dz;
   PetscInt     ind[24];
   PetscScalar  weight[24];
   PetscScalar3 vel;



 
   ierr = VecGetOwnershipRange(*x0,&Istart,&Iend);CHKERRQ(ierr); 
   Lsize = Iend - Istart; 

   PetscScalar xlast,ylast,zlast; 
   PetscScalar *xarray, *yarray,*zarray;  
   PetscScalar steplen;

   solsize   = cdinfo->nofgrid->x + cdinfo->nofgrid->y + cdinfo->nofgrid->z + 24;
   endx      = cdinfo->r->x;
   // add this to prevent the last grid problem..
   endx      = endx - cdinfo->gsize->x;
  
   for(i=0;i<24;i++){ VecSetValue(v,solsize-i,0, INSERT_VALUES);}
   ierr = VecGetArray(*x0,&xarray);CHKERRQ(ierr);
   ierr = VecGetArray(*y0,&yarray);CHKERRQ(ierr);
   ierr = VecGetArray(*z0,&zarray);CHKERRQ(ierr);

 
   ierr = VecCreateMPI(PETSC_COMM_WORLD,Lsize,PETSC_DETERMINE,&dydx);CHKERRQ(ierr); 
   ierr = VecCreateMPI(PETSC_COMM_WORLD,Lsize,PETSC_DETERMINE,&dzdx);CHKERRQ(ierr); 

   ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,solsize,&Ivec1);CHKERRQ(ierr); 
   ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,solsize,&Ivec2);CHKERRQ(ierr);
   ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,solsize,&Ivec3);CHKERRQ(ierr);
    
   ierr = VecGetOwnershipRange(Ivec1,&Ivstart,&Ivend);CHKERRQ(ierr);
   for(i=0;i<cdinfo->nofgrid->x;i++){
         if(i>=Ivstart&&i<Ivend){ierr = VecSetValue(Ivec1,i,1, INSERT_VALUES);}
                   } 
   for(i=cdinfo->nofgrid->x;i<cdinfo->nofgrid->x+cdinfo->nofgrid->y;i++){
         if(i>=Ivstart&&i<Ivend){ierr = VecSetValue(Ivec2,i,1, INSERT_VALUES);}
                   } 
   for(i=cdinfo->nofgrid->x+cdinfo->nofgrid->y;i<solsize-24;i++) {
         if(i>=Ivstart&&i<Ivend){ierr = VecSetValue(Ivec3,i,1, INSERT_VALUES);}
                   } 

   ierr = VecAssemblyBegin(Ivec1);
   ierr = VecAssemblyEnd(Ivec1);
   ierr = VecAssemblyBegin(Ivec2);
   ierr = VecAssemblyEnd(Ivec2);
   ierr = VecAssemblyBegin(Ivec3);
   ierr = VecAssemblyEnd(Ivec3);
  
// Now let's begin the streamline    
   ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,Lsize,solsize,(PetscInt)(40*cdinfo->gridspam->u.x),PETSC_NULL,&rowmat);CHKERRQ(ierr);
   for(j=0;j<Lsize;j++){
       xc = *(xarray+j);
       yc = *(yarray+j);
       zc = *(zarray+j);
       xlast   = xc;
       ylast   = yc;
       zlast   = zc;
       i       = 0;

       // test boundary
       //if(pow(yc-cdinfo->l->y,2)<1e-4) {xc=endx;}
       //if(pow(yc-cdinfo->r->y,2)<1e-4) {xc=endx;}
       //if(pow(zc-cdinfo->l->z,2)<1e-4) {xc=endx;}
       //if(pow(zc-cdinfo->l->z,2)<1e-4) {xc=endx;}    
    
       while(xc<endx && i<n-1){      
             RK2(&xlast, &ylast, &zlast, &v, cdinfo, &rowmat, j,dt);  
   
             steplen = sqrt(pow(xc-xlast,2) + pow(yc-ylast,2)+pow(zc-zlast,2));
        
             if(steplen >=2*(endx)/n||trycount>trycountmax){
                    xc = xlast;
                    yc = ylast;
                    zc = zlast;
                    trycount = 0;
                    i++;
             }else{
                    trycount++; 
             }

      }

 
     
     veval(xc, yc, zc, &v,cdinfo, &vel, ind, weight);
  
  
     VecSetValue(dydx,Istart+j,(PetscScalar)vel.y/vel.x, INSERT_VALUES); 
     VecSetValue(dzdx,Istart+j,(PetscScalar)vel.z/vel.x, INSERT_VALUES); 
     // Replace the locations by the end points 
     VecSetValue(*x0,Istart+j,xc, INSERT_VALUES); 
     VecSetValue(*y0,Istart+j,yc, INSERT_VALUES); 
     VecSetValue(*z0,Istart+j,zc, INSERT_VALUES); 
//fprintf(stdout,"%d %d pj= %d\n ",Istart,Iend,Istart+j);

} 

      //  x0,y0,z0 are sequential matrices!

      ierr = MatAssemblyBegin(rowmat,MAT_FINAL_ASSEMBLY);
      ierr = MatAssemblyEnd(rowmat,MAT_FINAL_ASSEMBLY);

  
      MatMerge(PETSC_COMM_WORLD,rowmat,PETSC_DECIDE,MAT_INITIAL_MATRIX,&weightmat);

     // MatDuplicate(weightmat,MAT_COPY_VALUES,&Dx); 
      MatDuplicate(weightmat,MAT_COPY_VALUES,&Dy);
      MatDuplicate(weightmat,MAT_COPY_VALUES,&Dz);
 
 
      // Dy = Dy + dy/dx * Dx
      //MatDiagonalScale(Dx,dydx,Ivec1);
      MatDiagonalScale(Dy,PETSC_NULL,Ivec2); 
     // MatAXPY(Dy,1,Dx,SAME_NONZERO_PATTERN);
      // Dz = Dz + dz/dx * Dx
      //MatDuplicate(weightmat,MAT_COPY_VALUES,&Dx); 
      //MatDiagonalScale(Dx,dzdx,Ivec1);
      MatDiagonalScale(Dz,PETSC_NULL,Ivec3);  
     // MatAXPY(Dz,1,Dx, SAME_NONZERO_PATTERN);




      *Dydv = Dy;
      *Dzdv = Dz;

   ierr = VecDestroy(dydx);CHKERRQ(ierr);
   ierr = VecDestroy(dzdx);CHKERRQ(ierr);
   ierr = VecDestroy(Ivec1);CHKERRQ(ierr);
   ierr = VecDestroy(Ivec2);CHKERRQ(ierr);
   ierr = VecDestroy(Ivec3);CHKERRQ(ierr);
   ierr = MatDestroy(weightmat);CHKERRQ(ierr);
   //ierr = MatDestroy(Dx);CHKERRQ(ierr);




return 0;
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
int StreamlineMaprr(Vec *x0, Vec *y0, Vec *z0, Vec v, CoordInfo *cdinfo,PetscScalar dt, PetscInt n ,PetscInt trycountmax) 
{
   PetscErrorCode ierr;
   PetscScalar    endx,xc,yc,zc; 
   PetscInt       solsize ,trycount = 0;
   PetscInt       i,j,Istart,Iend,Lsize;
   Mat            rowmat;
   PetscInt     ind[24];
   PetscScalar  weight[24];
   PetscScalar3 vel;
   
 
   ierr = VecGetOwnershipRange(*x0,&Istart,&Iend);CHKERRQ(ierr); 
   Lsize = Iend - Istart; 

   PetscScalar xlast,ylast,zlast; 
   PetscScalar *xarray, *yarray,*zarray;  
   PetscScalar steplen;

   solsize   = cdinfo->nofgrid->x + cdinfo->nofgrid->y + cdinfo->nofgrid->z + 24;
   endx      = cdinfo->r->x;

   for(i=0;i<24;i++){ VecSetValue(v,solsize-i,0, INSERT_VALUES);}

   // add this to prevent the last grid problem..
   endx      = endx - cdinfo->gsize->x;
  
   ierr = VecGetArray(*x0,&xarray);CHKERRQ(ierr);
   ierr = VecGetArray(*y0,&yarray);CHKERRQ(ierr);
   ierr = VecGetArray(*z0,&zarray);CHKERRQ(ierr);

 



// Now let's begin the streamline    
   
  // ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,Lsize,solsize,(PetscInt)(solsize*1.0/100),PETSC_NULL,&rowmat);CHKERRQ(ierr);

   for(j=0;j<Lsize;j++){
       xc = *(xarray+j);
       yc = *(yarray+j);
       zc = *(zarray+j);
       xlast   = xc;
       ylast   = yc;
       zlast   = zc;
       i       = 0;

       // test boundary
       //if(pow(yc-cdinfo->l->y,2)<1e-4) {xc=endx;}
       //if(pow(yc-cdinfo->r->y,2)<1e-4) {xc=endx;}
       //if(pow(zc-cdinfo->l->z,2)<1e-4) {xc=endx;}
       //if(pow(zc-cdinfo->l->z,2)<1e-4) {xc=endx;}    


    
       while(xc<endx && i<n-1){  
            
             RK2(&xlast, &ylast, &zlast, &v, cdinfo, PETSC_NULL, j,dt);    
             steplen = sqrt(pow(xc-xlast,2) + pow(yc-ylast,2)+pow(zc-zlast,2));
        
             if(steplen >=2*(endx)/n||trycount>trycountmax){
                    xc = xlast;
                    yc = ylast;
                    zc = zlast;
                    trycount = 0;
                    i++;
             }else{
                    trycount++; 
             }
         
         /*
             if(steplen >=2*(endx)/n){
                    xc = xlast;
                    yc = ylast;
                    zc = zlast;
                    trycount = 0;
                    i++;
             }else{
                 if(trycount>trycountmax){ 
                     xc = endx;
                     yc = *(yarray+j);
                     zc = *(zarray+j);
                    }

                    trycount++; 
             }
          */

      }

   
     
     veval(xc, yc, zc, &v,cdinfo, &vel, ind, weight);
  
 
     // Replace the locations by the end points 
     VecSetValue(*x0,Istart+j,xc, INSERT_VALUES); 
     VecSetValue(*y0,Istart+j,yc, INSERT_VALUES); 
     VecSetValue(*z0,Istart+j,zc, INSERT_VALUES); 

} 



  //   ierr = MatDestroy(rowmat);CHKERRQ(ierr);







return 0;
}



