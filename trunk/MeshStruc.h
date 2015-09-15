typedef struct
{
   PetscScalar x;
   PetscScalar y;
   PetscScalar z; 
} PetscScalar3;

typedef struct
{
   PetscInt x;
   PetscInt y;
   PetscInt z; 
} PetscInt3;

typedef struct
{
   PetscInt x;
   PetscInt y;
   PetscInt z;
   PetscInt w; 
} PetscInt4;

typedef struct
{
  PetscInt3 u;
  PetscInt3 v;
  PetscInt3 w;
} PetscInt33;



typedef struct 
{
  PetscInt     *dim;
  PetscInt3    *l;
  PetscInt3    *r;
  PetscInt3    *grid;
  PetscScalar3 *gsize;
  PetscInt3    *period;
 
  
  PetscInt33   *gridspam;
  PetscInt4    *nofgrid;

} CoordInfo;


typedef struct
{
 Mat Mat1;
 Mat Mat2;
 Mat Mat3; 
 Mat Mat4;
 Mat Mat5;
 Mat Mat6; 

} Mat6;

typedef struct
{
 Vec vec1;
 Vec vec2;
 Vec vec3; 
 Vec vec4;
 Vec vec5;
 Vec vec6;

} Vec6;

