//
// Created by Hell on 21.05.2018.
//

#include "suggested_work2.h"
#include "../legacy/myMatrices.h"

void getKrylov(StorageFormat *fmt, double *A, double **V, double **H, int j) {
    fmt->matvec(fmt, A, V[j], V[j + 1]);
    double *hj = H[j];
    double *w = V[j + 1];
    size_t n = (size_t) fmt->nx;
    for (int i = 0; i <= j; ++i) {
        hj[i] = dot(V[i], w, n);
        for (int k = 0; k < n; ++k) {
            w[k] -= hj[i] * V[i][k];
        }
    }
    hj[j + 1] = getNorm(w, n);
    constMult(w, 1 / hj[j + 1], w, n);
}

void applyGivensRotations(int j, double *c, double *s, double **H, double *g) {
    for (int k = 1; k < j + 1; ++k) {
        double t = c[k - 1] * H[j][k - 1] + s[k - 1] * H[j][k];
        H[j][k] = -s[k - 1] * H[j][k - 1] + c[k - 1] * H[j][k];
        H[j][k - 1] = t;
    }
    double sq = sqrt(H[j][j] * H[j][j] + H[j][j + 1] * H[j][j + 1]);
    c[j] = H[j][j] / sq;
    s[j] = H[j][j + 1] / sq;
    H[j][j] = sq;
//        H[j][j] = c[j]*H[j][j]+s[j]*H[j][j+1];
    g[j + 1] = -s[j] * g[j];
    g[j] = c[j] * g[j];
}

void backwardSubstitution(double **A, double *b, double *x, int size) {
    double sum = 0;
    x[size - 1] = b[size - 1] / A[size - 1][size - 1];
    for (int i = size - 2; i >= 0; --i) {
        for (int j = i + 1; j < size; ++j) {
            sum += A[j][i] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
        sum = 0;
    }
}

void assembleSolution(double **V, double *y, double *x0, double *xm, int n, int m) {
    //Assembling solution
    for (int i = 0; i < n; ++i) {
        xm[i] = x0[i];
        for (int j = 0; j < m; ++j) {
            xm[i] += y[j] * V[j][i];
        }
    }
}

void GMRES_malloc(int n, int m, double **r, double ***H, double ***V, double **g, double **y, double **c, double **s,
                  double **xm) {
    //error
    *g = newArray(m + 1);
    *y = newArray(m);
    *r = newArray(n);
    *xm = newArray(n);
    //Givens rotations
    *c = newArray(m);
    *s = newArray(m);
    //EigenMatrix and Hessenberg matrix
    *V = newMatrix(n, m + 1);
    *H = newMatrix(m + 1, m);
}

void GMRES_free(double *r, double **H, double **V, double *g, double *y, double *c, double *s, double *xm) {
    free(g);
    free(y);
    free(c);
    free(s);
    free(r);
    free(xm);
    freeMatrix(V);
    freeMatrix(H);
}

void GMRES_init(StorageFormat *fmt, double *A, double *b, double *x, double *r, double **V, double *g) {
    int n = fmt->nx;
    fmt->matvec(fmt, A, x, r);
    subtr(b, r, r, n);
    *g = getNorm(r, n);
    constMult(r, 1 / (*g), *V, n);
}

//GMRES-body
void GMRES_core(StorageFormat *fmt, double *A, double *x0, double **V, double **H,
                double *y, double *g, double *c, double *s, int m, double *xm, double *rsdl) {
    for (int j = 0; j < m; j++) {
        getKrylov(fmt, A, V, H, j);
        applyGivensRotations(j, c, s, H, g);
    }
    //Residual
    *rsdl = fabs(g[m]);
    backwardSubstitution(H, g, y, m);
    assembleSolution(V, y, x0, xm, fmt->nx, m);
}

//interface methods
double *GMRES_full(StorageFormat *fmt, double *A, double *x0, double *b, int m) {
    int n = fmt->nx;
    double *g, *y, *c, *s, **V, **H, *r, *xm, ro;
    double *x = newDoubleArray(n);
    GMRES_malloc(n, m, &r, &H, &V, &g, &y, &c, &s, &xm);
    GMRES_init(fmt, A, b, x0, r, V, g);
    GMRES_core(fmt, A, x0, V, H, y, g, c, s, m, xm, &ro);
    printf("Rsdl=%lg\n", ro);
    memcpy(x, xm, n * sizeof(double));
    GMRES_free(r, H, V, g, y, c, s, xm);
    return x;
}

double *GMRES_restarted(StorageFormat *fmt, double *A, double *x0, double *b, int m, double tol) {
    int n = fmt->nx;
    double *g, *y, *c, *s, **V, **H, *r, *xm, ro;
    double *x = newDoubleArray(n);
    memcpy(x, x0, n * sizeof(double));

    GMRES_malloc(n, m, &r, &H, &V, &g, &y, &c, &s, &xm);
    GMRES_init(fmt, A, b, x0, r, V, g);

    ro = g[0];
    while (ro > tol) {
        GMRES_core(fmt, A, x, V, H, y, g, c, s, m, xm, &ro);
        double *tmp = x;
        x = xm;
        xm = tmp;
        printf("Rsdl=%lg\n", ro);
        GMRES_init(fmt, A, b, x, r, V, g);
    }
    GMRES_free(r, H, V, g, y, c, s, xm);
    return x;
}

//CG METHOD

double * CG_body(StorageFormat *fmt, double *A, double *x0, double *b, double tol) {
    int n = fmt->nx;
    double *x = newDoubleArray(n);
    memcpy(x, x0, n * sizeof(double));
    double *r_m = newDoubleArray(n);
    double *p_m = newDoubleArray(n);
    double *Ap_m = newDoubleArray(n);
    fmt->matvec(fmt, A, x0, r_m);
    subtr(b, r_m, p_m, n);
    double r_m_norm_sq = dot(r_m, r_m, n);
    double ro = 1;
    double alpha;
    tol *= tol;
    while (ro > tol) {
        fmt->matvec(fmt, A, p_m, Ap_m); //Apm = A*pm
        alpha = r_m_norm_sq / dot(Ap_m, p_m, n);
        for (int i = 0; i < n; ++i) {
            x[i] = x[i] + alpha * p_m[i];
            r_m[i] = r_m[i] - alpha * Ap_m[i];
        }
        ro = dot(r_m, r_m, n);
        alpha = ro / r_m_norm_sq;
        for (int i = 0; i < n; ++i) {
            p_m[i] = r_m[i] + alpha * p_m[i];
        }
        r_m_norm_sq = ro;
    }
    return x;
}