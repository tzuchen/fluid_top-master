/*
Header file for functions in PetscReceive.c
*/

int VecReceive(PetscViewer, VecType, Vec*);
int MatReceiveTranspose(PetscViewer, MatType, Mat*);
int ScalarReceive(PetscViewer, PetscScalar**);
int IntReceive(PetscViewer, PetscInt**);

// function for my own use

int CdinfoReceive(PetscViewer, CoordInfo**);

int SparseView(Mat, PetscViewer);
int SparseSave(Mat,char[]);
int VecSave(Vec, char[]);
int VecSavebin(Vec vec, char[]);

int MatGetDataLocal(Mat, PetscInt**, PetscInt**, PetscScalar**, PetscInt*);

// function to load 
int VecLoadfromFile(MPI_Comm,char[], VecType,Vec*);
int MatLoadTransposefromFile(MPI_Comm,char[], MatType, Mat*);
int CdinfoLoad(MPI_Comm,char[], CoordInfo**);
int IntLoadfromFile(MPI_Comm, char[], PetscInt**);
int ScalarLoadfromFile(MPI_Comm, char[], PetscScalar**);
