

int LargeVecCreate(Vec*,PetscInt,Vec*);
int LargeVecGetOwnershipRange(Vec*,PetscInt,PetscInt*, PetscInt*);
int LargeVecGetColumnOwnershipRange(Vec*,PetscInt,PetscInt,PetscInt*, PetscInt* );
int ISCreateGeneralWithIJ(MPI_Comm,Vec,Vec*,PetscInt, PetscInt,PetscInt , PetscInt*, PetscInt*,IS*, IS*);
int LargeVecScatterCreate(Vec* ,IS* ,Vec ,IS* ,VecScatter* ,PetscInt);
int LargeVecScatterBeginEnd(Vec*,Vec,InsertMode,ScatterMode,VecScatter*,PetscInt);
int ISArrayDestroy(IS[],PetscInt);
int VecScatterArrayDestroy(VecScatter[],PetscInt);
int VecArrayDestroy(Vec[],PetscInt);
