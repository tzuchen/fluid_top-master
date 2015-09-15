/* 
Header file for functions in FFTRoutine.c
*/


int FFT2D(Vec, Vec, PetscInt, PetscInt, PetscInt, PetscInt, PetscInt, PetscTruth, PetscScalar,char*);
int SkewConvert(Vec, Vec,PetscInt,PetscInt,PetscInt);
int FFTColumn(Vec, Vec, PetscInt, PetscInt);
int IFFTColumn(Vec, Vec, PetscInt, PetscInt);
int FFTRow(Vec, Vec, PetscInt, PetscInt, PetscInt, PetscInt, PetscScalar**, PetscInt, PetscScalar);
int FFTRow2(Vec, Vec, PetscInt, PetscInt, PetscInt, PetscInt, PetscScalar**, PetscInt, PetscScalar);
int FFTRow3(Vec, Vec, PetscInt, PetscInt, PetscInt, PetscInt, PetscScalar**, PetscInt, PetscScalar);
int FFTRow4(Vec, Vec, PetscInt, PetscInt, PetscInt, PetscInt, PetscScalar**, PetscInt, PetscScalar);
int WaveNumber(PetscInt i,PetscInt j,PetscInt n);
