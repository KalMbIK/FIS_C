//
// Created by Hell on 21.05.2018.
//

#include "suggested_work2.h"
#include "../legacy/myMatrices.h"

void getKrylov(StorageFormat *fmt, double *A, double **V, double **H, int j) {
    fmt->matvec(fmt,A,V[j],V[j+1]);
    double *hj = H[j];
    double *w = V[j+1];
    size_t n = (size_t)fmt->nx;
    for (int i = 0; i <= j; ++i) {
        hj[i] = dot(V[i],w,n);
        for (int k = 0; k < n; ++k) {
            w[k] -= hj[i]*V[i][k];
        }
    }
    hj[j+1] = getNorm(w,n);
    constMult(w,1/hj[j+1],w,n);
}

void backwardSubstitution(double **A, double *b, double *x, int size){
    double sum = 0;
    x[size-1] = b[size-1]/A[size-1][size-1];
    for (int i = size-2; i >= 0; --i) {
        for (int j = i+1; j < size; ++j) {
            sum+=A[j][i]*x[j];
        }
        x[i] = (b[i] - sum)/A[i][i];
        sum = 0;
    }
}

void GMRES(StorageFormat *fmt, double *A, double *x0, double **V, double **H,
           double *y, double *g, double *c, double *s, int m, double *xm, double *rsdl) {
    int n = fmt->nx;
    //GMRES-body (Gram-Schmidt + Givens rotations on the fly)
    for (int j = 0; j < m; j++){
        getKrylov(fmt,A,V,H,j);
        for (int k = 1; k < j; ++k) {
            H[j][k-1] = c[k-1]*H[j][k-1]+s[k-1]*H[j][k];
            H[j][k] = -s[k-1]*H[j][k-1]+c[k-1]*H[j][k];
        }
        double sq = sqrt(H[j][j]*H[j][j]+H[j][j+1]*H[j][j+1]);
        c[j] = H[j][j]/sq;
        s[j] = H[j][j+1]/sq;
        H[j][j] = sq;
//        H[j][j] = c[j]*H[j][j]+s[j]*H[j][j+1];
        g[j+1] = -s[j]*g[j];
        g[j] = c[j]*g[j];
    }
    //Residual
    *rsdl = fabs(g[m]);

    //Assembling solution
    backwardSubstitution(H,g,y,m);
//    double sum = 0;
//    for (int l = 0; l < m; ++l) {
//        double dotSum = 0;
//        for (int i = l; i < m; ++i) {
//            dotSum += H[i][l]*y[i];
//        }
//        sum += g[l]-dotSum;
//    }
//    printf("sum=%lg\n",sum);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            xm[i]=x0[i]+y[j]*V[j][i];
        }
    }
}

double * restartedGMRES(StorageFormat *fmt, double *A, double *b, double *x0, int m, double tol) {
    int n = fmt->nx;
    //_________________________
    //error
    double *g = newArray(m+1);
    double *y = newArray(m);
    //Givens rotations
    double *c = newArray(m);
    double *s = newArray(m);
    //EigenMatrix and Hessenberg matrix
    double **V = newMatrix(n,m+1);
    double **H = newMatrix(m+1,m);
    //_________________________
    double *x = newArray(n);
    double *r = newArray(n);
    double *xm = newArray(n);
    fmt->matvec(fmt,A,x0,r);
    subtr(b,r,r,n);
    double beta, ro = getNorm(r,n);
    constMult(r,1/ro,V[0],n);
    memcpy(x,x0,n);
    //beta = ||r0||
    g[0] = ro;
    while (ro > tol){
        GMRES(fmt, A, x, V, H, y, g, c, s, m, xm, &ro);
        double *tmp = x;
        x = xm;
        xm = tmp;
        printf("Rsdl=%lg\n",ro);
        fmt->matvec(fmt,A,x,r);
        subtr(b,r,r,n);
        beta = getNorm(r,n);
        constMult(r,1/beta,V[0],n);
//        printf("%lg\n",beta);
//        printf("%lg\n",getNorm(V[0],n));
        g[0] = beta;
    }
    free(g);
    free(y);
    free(c);
    free(s);
    free(r);
    free(xm);
    freeMatrix(V);
    freeMatrix(H);
    return x;
}