/*
Header file for iteration scheme

Tzu-Chen Liang  11-6-2006

*/

int BackwardAverage(Vec*, Vec*, PetscInt*, PetscScalar* , PetscInt, PetscInt , PetscInt, PetscInt, PetscInt, PetscScalar);
int BackwardAverageRL(Vec*, Vec*, PetscInt*, PetscScalar* , PetscInt, PetscInt , PetscInt, PetscInt, PetscInt, PetscScalar);
int BackwardCenter(Vec*, Vec*, PetscInt*, PetscScalar* , PetscInt, PetscInt , PetscInt, PetscInt, PetscInt, PetscScalar);
int Smoothing(Vec*,Vec*,PetscScalar*,PetscInt*,VecScatter*,PetscInt,DA,PetscInt,PetscInt);
int SmoothingRL(Vec*,Vec*,PetscScalar*,PetscInt*,VecScatter*,PetscInt,PetscInt,PetscInt);
int InverseStandardMap(PetscScalar*,PetscScalar*,PetscScalar);
int StandardMap(PetscScalar*,PetscScalar*,PetscScalar);
int InverseModifiedArnoldsCatMap(PetscScalar*,PetscScalar*);
int SkewSymmetricPoint(PetscInt*, PetscInt*, PetscInt);

int SkewSymmetricScatter(Vec*, PetscScalar*, PetscInt*,PetscInt,PetscInt, PetscInt, VecScatter*);
