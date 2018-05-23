//
// Created by Hell on 21.05.2018.
//

#ifndef FIS_C_SUGGESTED_WORK2_H
#define FIS_C_SUGGESTED_WORK2_H

#include "suggested_work1.h"

void getKrylov(StorageFormat *fmt, double *A, double **V, double **H, int j);
void backwardSubstitution(double **A, double *b, double *x, int size);
void GMRES(StorageFormat *fmt, double *A, double *x0, double **V, double **H, double *y, double *g, double *c, double *s,
           int m, double *xm, double *rsdl);
double * restartedGMRES(StorageFormat *fmt, double *A, double *b, double *x0, int m, double tol);

#endif //FIS_C_SUGGESTED_WORK2_H
