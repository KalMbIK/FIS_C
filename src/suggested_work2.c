//
// Created by Hell on 21.05.2018.
//

#include "suggested_work2.h"
#include "../legacy/myMatrices.h"

void getKrylov(StorageFormat *fmt, double *A, double **V, double **H, int j) {
//    fmt->matvec(fmt, A, V[j], V[j + 1]);
    fmt->matvec(fmt, A, V[j], fmt->Z);
    sparseDiag(fmt, A, fmt->Z, V[j + 1]);
//    sparseForwardSubstitution(fmt, A, fmt->Z, V[j + 1]);
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
//    *g = getNorm(r, n);
//    constMult(r, 1 / (*g), *V, n);
    sparseDiag(fmt, A, r, *V);
//    sparseForwardSubstitution(fmt, A, b, *V);
    *g = getNorm(*V, n);
    constMult(*V, 1 / (*g), *V, n);
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

void CG_malloc(int n, double **r_m, double **p_m, double **Ap_m) {
    double *x = newDoubleArray(3 * n);
    *r_m = x;
    *p_m = x + n;
    *Ap_m = x + n + n;
}

void CG_free(double *r_m, double *p_m, double *Ap_m) {
    free(r_m);
}

double *CG_body(StorageFormat *fmt, double *A, double *x0, double *b, double tol) {
    int n = fmt->nx;
    double *x = newDoubleArray(n);
    memcpy(x, x0, n * sizeof(double));

    int iters = 0;
    double qtol = tol * tol;
    double ro = 1;
    double alpha;
    double *r_m, *p_m, *Ap_m;

    CG_malloc(n, &r_m, &p_m, &Ap_m);
    //r0=b-Ax0, p=r0, (r0,r0)
    fmt->matvec(fmt, A, x0, r_m);
    subtr(b, r_m, r_m, n);
    memcpy(p_m, r_m, n * sizeof(double));
//    printDoubleArray(p_m,n);
    double r_m_norm_sq = dot(r_m, r_m, n);
    printf("%lg\n", r_m_norm_sq);
    double r_0_norm_sq = r_m_norm_sq, rel_error = 1;

    while (rel_error > qtol) {
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
        rel_error = r_m_norm_sq/r_0_norm_sq;
        if (iters % 1000 == 0){
            printf("Iter: %d, Rel_error=%lg\n", iters, rel_error);
        }
        iters++;
    }
    printf("---\nFINISHED\nIter: %d, Rel_error=%lg\n---\n", iters, rel_error);

    CG_free(r_m, p_m, Ap_m);
    return x;
}

//PRECONDITIONING
void sparseForwardSubstitution(StorageFormat *fmt, double *A, double *b, double *x) {
    x[0] = b[0] / A[0];
    int n = fmt->nx;
    for (int i = 1; i < n; ++i) {
        int i1 = fmt->I[i];
        int i2 = fmt->I[i + 1];
        double sum = 0;
        int j;
        for (j = i1; j < i2; ++j) {
            if (fmt->J[j] < i) {
                sum += A[j] * x[fmt->J[j]];
            } else {
                break;
            }
        }
        x[i] = (b[i] - sum) / A[j];
    }
}

void sparseDiag(StorageFormat *fmt, double *A, double *b, double *x) {
    int n = fmt->nx;
    for (int i = 0; i < n; ++i) {
        int i1 = fmt->I[i];
        int i2 = fmt->I[i + 1];
        for (int j = i1; j < i2; ++j) {
            if (fmt->J[j] == i) {
                x[i] = b[i] / A[j];
                break;
            }
        }
    }
}