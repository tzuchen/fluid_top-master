/*
Here are the functions to calcuate Markov matrix 
from a given flow velocity field

Tzu-Chen Liang 2-28-2007
*/

int MarkovMap(Vec, PetscInt, PetscInt, CoordInfo*,Mat*, Mat6*, Vec6*, PetscScalar**, PetscInt*);
int dMarkovMap(Vec*, Vec6, Mat6);

int MarkovMap2(Vec, PetscInt, PetscInt, CoordInfo*,Mat*);
int MarkovMap3(Vec, PetscInt, PetscInt, CoordInfo*,Mat*);
int MarkovMap4(Vec, PetscInt, PetscInt, CoordInfo*,Mat*);

int MarkovMap2d(Vec, PetscInt, PetscInt, CoordInfo*,Mat*,Mat*,Mat*);
int MarkovMap2dv(Vec, PetscInt, PetscInt, CoordInfo*,Mat*,Mat*,Mat*);
int dMarkovMap2d(Vec, PetscInt, PetscInt, CoordInfo*,Mat*, Mat*);
int dMarkovMap2df(Vec, PetscInt, PetscInt, CoordInfo*,Vec*);
