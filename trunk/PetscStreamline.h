/*
Headers of functions in PetscStreamline.c
All functions here are called by RK2 or 4. 

Tzu-Chen Liang   2-8-2006
*/


int  loc2ijk(PetscInt, PetscInt, PetscInt, PetscScalar, PetscScalar, PetscInt*, PetscScalar*, PetscInt);
int gloc2ijk(PetscInt3*, PetscScalar3*, PetscScalar, PetscScalar, PetscScalar, PetscInt*, PetscScalar*, PetscInt33*);
int  ijk2ind(PetscInt*, PetscInt3*, PetscInt3*, PetscInt*);
int gijk2ind(PetscInt*, PetscInt33*,PetscInt3*, PetscInt*);
int  cornerfind(PetscInt*, PetscInt*);
int gcornerfind(PetscInt*, PetscInt*);
int  lind2gind(PetscInt*, PetscInt, PetscInt, PetscInt);
int glind2gind(PetscInt*, PetscInt4*, PetscInt);
int  lininterp3d(PetscScalar*, PetscScalar*, PetscScalar*, PetscScalar*);
int glininterp3d(PetscScalar*, PetscScalar*, PetscScalar3*, PetscScalar*);
int veval(PetscScalar, PetscScalar, PetscScalar, Vec* ,CoordInfo*, PetscScalar3*,PetscInt* , PetscScalar*);
int RK2(PetscScalar*, PetscScalar*, PetscScalar*, Vec*, CoordInfo*,Mat* ,PetscInt ,PetscScalar);
int RK4(PetscScalar*, PetscScalar*, PetscScalar*, Vec*, CoordInfo*,Mat* ,PetscInt ,PetscScalar);
int StreamlineMap(Vec, Vec, Vec, Vec, CoordInfo*, Mat*, Mat*,PetscScalar, PetscInt, PetscInt,PetscScalar**, PetscScalar**, PetscScalar**);
int StreamlineMapr(Vec*, Vec*, Vec*, Vec, CoordInfo*, Mat*, Mat*,PetscScalar, PetscInt, PetscInt);

int StreamlineMaprr(Vec*, Vec*, Vec*, Vec, CoordInfo*,PetscScalar, PetscInt ,PetscInt );


