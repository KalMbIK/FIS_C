//
// Created by Hell on 21.05.2018.
//

#ifndef FIS_C_SUGGESTED_WORK2_H
#define FIS_C_SUGGESTED_WORK2_H

#include "suggested_work1.h"

struct SparseSolver{
    StorageFormat *fmt;
    double *A, *b, *x;

    double (*solve)(struct StorageFormat*, ...);
};
typedef struct SparseSolver SparseSolver;

void getKrylov(StorageFormat *fmt, double *A, double **V, double **H, int j);
void applyGivensRotations(int j, double* c, double* s, double** H, double* g);
void assembleSolution(double** V, double* y, double* x0, double* xm, int n, int m);
void backwardSubstitution(double **A, double *b, double *x, int size);

//GMRES_METHOD
void GMRES_malloc(int n, int m,
                  double **r, double ***H, double ***V, double **g, double **y, double **c, double **s, double **xm);
void GMRES_free(double *r, double **H, double **V, double *g, double *y, double *c, double *s, double *xm);

void GMRES_init(StorageFormat *fmt, double *A, double *b, double *x, double *r, double **V, double *g);

void GMRES_core(StorageFormat *fmt, double *A, double *x0,
                double **V, double **H, double *y, double *g, double *c, double *s, int m, double *xm, double *rsdl);
double *GMRES_full(StorageFormat *fmt, double *A, double *x0, double *b, int m);

double *GMRES_restarted(StorageFormat *fmt, double *A, double *x0, double *b, int m, double tol);

//CG METHOD
void CG_malloc(int n, double **r_m, double **p_m, double **Ap_m);

void CG_free(double *r_m, double *p_m, double *Ap_m);

double * CG_body(StorageFormat *fmt, double *A, double *x0, double *b, double tol);

//PRECONDITIONING
void sparseForwardSubstitution(StorageFormat *fmt, double *A, double *b, double *x);

void sparseDiag(StorageFormat *fmt, double *A, double* b, double* x);
#endif //FIS_C_SUGGESTED_WORK2_H
