static char help[] = "A Test for Socket programming.";

#include "petscksp.h"
#include "petscda.h"

#if defined(PETSC_NEEDS_UTYPE_TYPEDEFS)
/* Some systems have inconsistent include files that use but don't
   ensure that the following definitions are made */
typedef unsigned char   u_char;
typedef unsigned short  u_short;
typedef unsigned int    u_int;
typedef unsigned long   u_long;
#endif

#include <errno.h>
#if defined(PETSC_HAVE_STDLIB_H)
#include <stdlib.h>
#endif
#include <sys/types.h>
#include <ctype.h>
#if defined(PETSC_HAVE_MACHINE_ENDIAN_H)
#include <machine/endian.h>
#endif
#if defined(PETSC_HAVE_UNISTD_H)
#include <unistd.h>
#endif
#if defined(PETSC_HAVE_SYS_SOCKET_H)
#include <sys/socket.h>
#endif
#if defined(PETSC_HAVE_SYS_WAIT_H)
#include <sys/wait.h>
#endif
#if defined(PETSC_HAVE_NETINET_IN_H)
#include <netinet/in.h>
#endif
#if defined(PETSC_HAVE_NETDB_H)
#include <netdb.h>
#endif
#if defined(PETSC_HAVE_FCNTL_H)
#include <fcntl.h>
#endif
#if defined(PETSC_HAVE_STROPTS_H)
#include <stropts.h>
#endif
#if defined (PETSC_HAVE_IO_H)
#include <io.h>
#endif
#if defined(PETSC_HAVE_SYS_UTSNAME_H)
#include <sys/utsname.h>
#endif
#if defined(PETSC_HAVE_STRINGS_H)
#include <strings.h>
#endif
#if defined(PETSC_HAVE_STRING_H)
#include <string.h>
#endif
#if defined(PETSC_HAVE_WINSOCK2_H)
#include <Winsock2.h>
#endif
#if defined(PETSC_HAVE_WS2TCPIP_H)
#include <Ws2tcpip.h>
#endif
#include "src/sys/src/viewer/impls/socket/socket.h"
#include "petscfix.h"


#undef __FUNCT__
#define __FUNCT__ "main"
EXTERN PetscErrorCode SOCKConnect_Private(int);

int main(int argc,char **args)
{
 PetscErrorCode ierr;
 PetscInt       rank,size,i,Istart,Iend,m,n,s;
 PetscScalar    q;
 Vec		x;
 int 		t,portnumber;
 

  PetscInitialize(&argc,&args,(char *)0,help);
 
  MPI_Comm_size(PETSC_COMM_WORLD,&size);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  
 	  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
      	  ierr = VecSetSizes(x,PETSC_DECIDE,10);CHKERRQ(ierr);        	  
          ierr = VecSetType(x, VECMPI);CHKERRQ(ierr); 
          ierr = VecGetOwnershipRange(x,&Istart,&Iend);CHKERRQ(ierr);

     
          // Use the default portnumber. Notice only one machine (rank 0) is connected to the socket. 
          portnumber=5050;
          if (rank ==0){t = SOCKConnect_Private(portnumber);}

// Receive the data (again, only rank 0 does this)
///////////////////////////////////////////////////////////////////////
// A Vector
   if (rank==0){
          PetscInt VecSize;
          PetscScalar *p1;
          
         
          // The first one received is the size of the vector          
          ierr = PetscBinaryRead(t,&m,1,PETSC_DOUBLE);CHKERRQ(ierr);
          fprintf(stdout,"m = %d\n",m);
          VecSize = (PetscInt)m;
          p1 = malloc(2*VecSize*sizeof(PETSC_DOUBLE)); // not sure why we need a 2*
          
     
          // The second one received is 0 (not sure waht it means) 
          ierr = PetscBinaryRead(t,&n,1,PETSC_DOUBLE);CHKERRQ(ierr);
          fprintf(stdout,"n = %d\n",n);     
    
	
          ierr = PetscBinaryRead(t,p1,VecSize,PETSC_DOUBLE);CHKERRQ(ierr);
          for (i =0; i<VecSize; i++) {ierr = VecSetValue(x,i,*(p1+i),INSERT_VALUES);CHKERRQ(ierr);}          

          free(p1);
        };  

          
   // Scatter the received vec to other processors.      
           ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
	   ierr = VecAssemblyEnd(x);CHKERRQ(ierr);
           VecView(x,0);
           ierr = VecDestroy(x);CHKERRQ(ierr);
//////////////////////////////////////////////////////////////////////////
// A Scalar
if(rank ==0){
 ierr = PetscBinaryRead(t,&q,1,PETSC_DOUBLE);CHKERRQ(ierr);
 fprintf(stdout,"q = %f\n",q);    
         }

//////////////////////////////////////////////////////////////////////////
// An Integer
if(rank ==0){
 ierr = PetscBinaryRead(t,&s,1,PETSC_INT);CHKERRQ(ierr);
 fprintf(stdout,"s = %i\n",s);    
         }

//////////////////////////////////////////////////////////////////////////



  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}

static int listenport;
/*-----------------------------------------------------------------*/
extern PetscErrorCode establish(u_short);
#undef __FUNCT__  
#define __FUNCT__ "SOCKConnect_Private"
PetscErrorCode SOCKConnect_Private(int portnumber)
{
struct sockaddr_in isa; 
#if defined(PETSC_HAVE_ACCEPT_SIZE_T)
  size_t             i;
#else
  int                i;
#endif
  int                t;

/* open port*/
  listenport = establish((u_short) portnumber);
  if (listenport == -1) {
       fprintf(stdout,"RECEIVE: unable to establish port\n");
       return -1;
  }

/* wait for someone to try to connect */
  i = sizeof(struct sockaddr_in);
  if ((t = accept(listenport,(struct sockaddr *)&isa,(socklen_t *)&i)) < 0) {
     fprintf(stdout,"RECEIVE: error from accept\n");
     return(-1);
  }
  close(listenport);  
  return(t);
}
/*-----------------------------------------------------------------*/

/*-----------------------------------------------------------------*/
#define MAXHOSTNAME 100
#undef __FUNCT__  
#define __FUNCT__ "establish"
PetscErrorCode establish(u_short portnum)
{
  char               myname[MAXHOSTNAME+1];
  int                s;
  PetscErrorCode     ierr;
  struct sockaddr_in sa;  
  struct hostent     *hp;
#if defined(PETSC_HAVE_UNAME)
  struct utsname     utname;
#elif defined(PETSC_HAVE_GETCOMPUTERNAME)
  int                namelen=MAXHOSTNAME;
#endif

  /* Note we do not use gethostname since that is not POSIX */
#if defined(PETSC_HAVE_GETCOMPUTERNAME)
  GetComputerName((LPTSTR)myname,(LPDWORD)&namelen);
#elif defined(PETSC_HAVE_UNAME)
  uname(&utname);
  strncpy(myname,utname.nodename,MAXHOSTNAME);
#endif
#if defined(PETSC_HAVE_BZERO)
  bzero(&sa,sizeof(struct sockaddr_in));
#else
  memset(&sa,0,sizeof(struct sockaddr_in));
#endif
  hp = gethostbyname(myname);


  if (!hp) {
     fprintf(stdout,"RECEIVE: error from gethostbyname\n");
     return(-1);
  }

  sa.sin_family = hp->h_addrtype; 
  sa.sin_port = htons(portnum); 
 
  if ((s = socket(AF_INET,SOCK_STREAM,0)) < 0) {
     fprintf(stdout,"RECEIVE: error from socket\n");
     return(-1);
  }
  {
  int optval = 1; /* Turn on the option */
  ierr = setsockopt(s,SOL_SOCKET,SO_REUSEADDR,(char *)&optval,sizeof(optval));
  }

  while (bind(s,(struct sockaddr*)&sa,sizeof(sa)) < 0) {
#if defined(PETSC_HAVE_WSAGETLASTERROR)
    ierr = WSAGetLastError();
    if (ierr != WSAEADDRINUSE) {
#else
    if (errno != EADDRINUSE) { 
#endif
      close(s);
      fprintf(stdout,"RECEIVE: error from bind\n");
      return(-1);
    }
    close(listenport); 
  }
//fprintf(stdout,"s = %i\n",s);
  listen(s,0);
  return(s);
}


