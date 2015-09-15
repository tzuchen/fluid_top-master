
#include <stdio.h>
#include "petscsys.h"
#include "src/sys/viewer/impls/socket/socket.h"
#include "mex.h"


/*
   TAKEN from src/sys/src/fileio/sysio.c The swap byte routines are 
  included here because the Matlab programs that use this do NOT
  link to the PETSc libraries.
*/
#include <errno.h>
#if defined(PETSC_HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if !defined(PETSC_WORDS_BIGENDIAN)
/*
  SYByteSwapInt - Swap bytes in an integer
*/
#undef __FUNCT__  
#define __FUNCT__ "SYByteSwapInt"
void SYByteSwapInt(int *buff,int n)
{
  int  i,j,tmp;
  char *ptr1,*ptr2 = (char*)&tmp;
  for (j=0; j<n; j++) {
    ptr1 = (char*)(buff + j);
    for (i=0; i<sizeof(int); i++) {
      ptr2[i] = ptr1[sizeof(int)-1-i];
    }
    buff[j] = tmp;
  }
}
/*
  SYByteSwapShort - Swap bytes in a short
*/
#undef __FUNCT__  
#define __FUNCT__ "SYByteSwapShort"
void SYByteSwapShort(short *buff,int n)
{
  int   i,j;
  short tmp;
  char  *ptr1,*ptr2 = (char*)&tmp;
  for (j=0; j<n; j++) {
    ptr1 = (char*)(buff + j);
    for (i=0; i<sizeof(short); i++) {
      ptr2[i] = ptr1[sizeof(int)-1-i];
    }
    buff[j] = tmp;
  }
}
/*
  SYByteSwapScalar - Swap bytes in a double
  Complex is dealt with as if array of double twice as long.
*/
#undef __FUNCT__  
#define __FUNCT__ "SYByteSwapScalar"
void SYByteSwapScalar(PetscScalar *buff,int n)
{
  int    i,j;
  double tmp,*buff1 = (double*)buff;
  char   *ptr1,*ptr2 = (char*)&tmp;
#if defined(PETSC_USE_COMPLEX)
  n *= 2;
#endif
  for (j=0; j<n; j++) {
    ptr1 = (char*)(buff1 + j);
    for (i=0; i<sizeof(double); i++) {
      ptr2[i] = ptr1[sizeof(double)-1-i];
    }
    buff1[j] = tmp;
  }
}
#endif

#undef __FUNCT__  
#define __FUNCT__ "PetscBinaryRead"
/*
    PetscBinaryRead - Reads from a binary file.

  Input Parameters:
.   fd - the file
.   n  - the number of items to read 
.   type - the type of items to read (PETSC_INT or PETSC_SCALAR)

  Output Parameters:
.   p - the buffer

  Notes: does byte swapping to work on all machines.
*/
PetscErrorCode PetscBinaryRead(int fd,void *p,int n,PetscDataType type)
{

  int  maxblock,wsize,err;
  char *pp = (char*)p;
#if !defined(PETSC_WORDS_BIGENDIAN)
  int  ntmp = n; 
  void *ptmp = p; 
#endif

  maxblock = 65536;
  if (type == PETSC_INT)         n *= sizeof(int);
  else if (type == PETSC_SCALAR) n *= sizeof(PetscScalar);
  else if (type == PETSC_SHORT)  n *= sizeof(short);
  else printf("PetscBinaryRead: Unknown type");
  
  while (n) {
    wsize = (n < maxblock) ? n : maxblock;
    err = read(fd,pp,wsize);
#if !defined(PETSC_MISSING_ERRNO_EINTR)
    if (err < 0 && errno == EINTR) continue;
#endif
    if (!err && wsize > 0) return 1;
    if (err < 0) {
      perror("error reading");
      return err;
    }
    n  -= err;
    pp += err;
  }
#if !defined(PETSC_WORDS_BIGENDIAN)
  if (type == PETSC_INT) SYByteSwapInt((int*)ptmp,ntmp);
  else if (type == PETSC_SCALAR) SYByteSwapScalar((PetscScalar*)ptmp,ntmp);
  else if (type == PETSC_SHORT) SYByteSwapShort((short*)ptmp,ntmp);
#endif

  return 0;
}


////////////////////////////////////////////////////////////////////////////

 PetscErrorCode PetscBinaryWrite(int fd,void *p,PetscInt n,PetscDataType type,PetscTruth istemp)
 {
   char           *pp = (char*)p;
   int            err,wsize;
   size_t         m = (size_t)n,maxblock=65536;
 #if !defined(PETSC_WORDS_BIGENDIAN) || (PETSC_SIZEOF_INT == 8) ||  defined(PETSC_USE_64BIT_INT)
   void           *ptmp = p;
 #endif

   //if (n < 0) SETERRQ1(PETSC_ERR_ARG_OUTOFRANGE,"Trying to write a negative amount of data %D",n);
   if (!n) return(0);

   if (type == PETSC_INT){
     m   *= sizeof(int);
 #if (PETSC_SIZEOF_INT == 8) || defined(PETSC_USE_64BIT_INT)
     PetscInt   *p_int = (PetscInt*)p,i;
     int *p_short;
     PetscMalloc(m,&pp);
     ptmp    = (void*)pp;
     p_short = (Int*)pp;

     for (i=0; i<n; i++) {
       p_short[i] = (int) p_int[i];
     }
     istemp = PETSC_TRUE;
 #endif
   }
   else if (type == PETSC_SCALAR)  m *= sizeof(PetscScalar);
   else if (type == PETSC_DOUBLE)  m *= sizeof(double);
   else if (type == PETSC_SHORT)   m *= sizeof(short);
   else if (type == PETSC_CHAR)    m *= sizeof(char);
   else if (type == PETSC_ENUM)    m *= sizeof(PetscEnum);
   else if (type == PETSC_TRUTH)   m *= sizeof(PetscTruth);
   //else if (type == PETSC_LOGICAL) m = PetscBTLength(m)*sizeof(char);
   else {}//SETERRQ(PETSC_ERR_ARG_OUTOFRANGE,"Unknown type");

 #if !defined(PETSC_WORDS_BIGENDIAN)
   if      (type == PETSC_INT)    {PetscByteSwapInt((int*)ptmp,n);}
   else if (type == PETSC_ENUM)   {PetscByteSwapInt((int*)ptmp,n);}
   else if (type == PETSC_TRUTH)  {PetscByteSwapInt((int*)ptmp,n);}
   else if (type == PETSC_SCALAR) {PetscByteSwapScalar((PetscScalar*)ptmp,n);}
   else if (type == PETSC_DOUBLE) {PetscByteSwapDouble((double*)ptmp,n);}
   else if (type == PETSC_SHORT)  {PetscByteSwapShort((short*)ptmp,n);}
 #endif

   while (m) {
     wsize = (m < maxblock) ? m : maxblock;
     err = write(fd,pp,wsize);
     if (err < 0 && errno == EINTR) continue;
     //if (err != wsize) SETERRQ(PETSC_ERR_FILE_WRITE,"Error writing to file.");
     m -= wsize;
     pp += wsize;
   }

 #if !defined(PETSC_WORDS_BIGENDIAN) && !(PETSC_SIZEOF_INT == 8) && !defined(PETSC_USE_64BIT_INT)
   if (!istemp) {
     if      (type == PETSC_SCALAR) {PetscByteSwapScalar((PetscScalar*)ptmp,n);}
     else if (type == PETSC_SHORT)  {PetscByteSwapShort((short*)ptmp,n);}
     else if (type == PETSC_INT)    {PetscByteSwapInt((int*)ptmp,n);}
     else if (type == PETSC_ENUM)   {PetscByteSwapInt((int*)ptmp,n);}
     else if (type == PETSC_TRUTH)  {PetscByteSwapInt((int*)ptmp,n);}
   }
 #endif

 #if (PETSC_SIZEOF_INT == 8) || defined(PETSC_USE_64BIT_INT)
   if (type == PETSC_INT){
     PetscFree(ptmp);
   }
 #endif
   return(0);
 }



/////////////////////////////////////////////////////////////////////////////////////////////////

  PetscErrorCode PETSC_DLLEXPORT PetscByteSwapInt(int *buff,PetscInt n)
 {
    PetscInt  i,j,tmp = 0;
    PetscInt  *tptr = &tmp;                /* Need to access tmp indirectly to get */
    char      *ptr1,*ptr2 = (char*)&tmp; /* arround the bug in DEC-ALPHA g++ */
  
    for (j=0; j<n; j++) {
      ptr1 = (char*)(buff + j);
      for (i=0; i<(int)sizeof(int); i++) {
        ptr2[i] = ptr1[sizeof(int)-1-i];
      }
      buff[j] = *tptr;
    }
    return(0);
  }
  /* --------------------------------------------------------- */
  /*
    PetscByteSwapShort - Swap bytes in a short
  */
  PetscErrorCode PETSC_DLLEXPORT PetscByteSwapShort(short *buff,PetscInt n)
  {
    PetscInt   i,j;
    short      tmp;
    short      *tptr = &tmp;           /* take care pf bug in DEC-ALPHA g++ */
    char       *ptr1,*ptr2 = (char*)&tmp;

    for (j=0; j<n; j++) {
      ptr1 = (char*)(buff + j);
      for (i=0; i<(int) sizeof(short); i++) {
        ptr2[i] = ptr1[sizeof(int)-1-i];
      }
      buff[j] = *tptr;
    }
    return(0);
  }
  /* --------------------------------------------------------- */
  /*
    PetscByteSwapScalar - Swap bytes in a double
    Complex is dealt with as if array of double twice as long.
  */
  PetscErrorCode PETSC_DLLEXPORT PetscByteSwapScalar(PetscScalar *buff,PetscInt n)
  {
    PetscInt  i,j;
    PetscReal tmp,*buff1 = (PetscReal*)buff;
    PetscReal *tptr = &tmp;          /* take care pf bug in DEC-ALPHA g++ */
    char      *ptr1,*ptr2 = (char*)&tmp;

  #if defined(PETSC_USE_COMPLEX)
    n *= 2;
  #endif
    for (j=0; j<n; j++) {
      ptr1 = (char*)(buff1 + j);
      for (i=0; i<(int) sizeof(PetscReal); i++) {
        ptr2[i] = ptr1[sizeof(PetscReal)-1-i];
      }
      buff1[j] = *tptr;
    }
    return(0);
 }
 /* --------------------------------------------------------- */
 /*
   PetscByteSwapDouble - Swap bytes in a double
 */
 PetscErrorCode PETSC_DLLEXPORT PetscByteSwapDouble(double *buff,PetscInt n)
 {
   PetscInt i,j;
   double   tmp,*buff1 = (double*)buff;
   double   *tptr = &tmp;          /* take care pf bug in DEC-ALPHA g++ */
   char     *ptr1,*ptr2 = (char*)&tmp;

   for (j=0; j<n; j++) {
     ptr1 = (char*)(buff1 + j);
     for (i=0; i<(int) sizeof(double); i++) {
       ptr2[i] = ptr1[sizeof(double)-1-i];
     }
     buff1[j] = *tptr;
   }
   return(0);
 }
 //#endif
 /* --------------------------------------------------------- */
 //@C


//////////////////////////////////////////////////////////////////////////////////

 PetscErrorCode  PetscBinaryOpen(const char name[],PetscFileMode mode,int *fd)
 {

   *fd = creat(name,"w",0666);

   return(0);
 }

//////////////////////////////////////////////////////////////////////////////////

 PetscErrorCode  PetscBinaryClose(int fd)
 {
   close(fd);
   return(0);
 }






