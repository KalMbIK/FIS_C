#include "suggested_work2.h"

void testMethods(int *argc, char ***argv);
double* readVec_from_vecfile(char* vecfile, int* size);

int main(int argc, char* argv[]) {
    testMethods(&argc, &argv);
    return 0;
}

void testMethods(int *argc, char ***argv){
    StorageFormat coo, csr;
    double *V = mm_read_mtx_to_COO((*argv)[3], &coo);
    formatConvertCOOtoCSR(&coo, &csr);
    csr.Z = newArray(csr.nx);
    size_t n = (size_t)csr.nx;

    double *x0 = newDoubleArray(n);
    double *b = newDoubleArray(n);
    double *x_star = ones(n);
    csr.matvec(&csr, V, x_star, b);

/*MATVEC TEST*/
//    int n2;
//    double *x = readVec_from_vecfile((*argv)[4], &n2);
//    if (n==n2){
//        subtr(b,b2,x_star,n);
//    }

/*TEST GMRES_full\GMRES_restarted\CG*/
//    double *x = GMRES_full(&csr, V, x0, b, atoi((*argv)[1]));
//    double  *x = GMRES_restarted(&csr, V, x0, b, atoi((*argv)[1]), atof((*argv)[2]));
    double  *x = CG_body(&csr, V, x0, b, atof((*argv)[2]));
    subtr(x,x_star,x_star,n);
    printDoubleArray(x, n);
    printf("||error||=%lg\n",getNorm(x_star,n));

    free(csr.Z);
    formatFree(&coo);
    formatFree(&csr);
    free(V);
    free(x0);
    free(x);
    free(x_star);
    free(b);
}

double* readVec_from_vecfile(char* vecfile, int* size){
    FILE* file = fopen(vecfile,"r");
    fscanf(file,"%d",size);
    double * x = newDoubleArray(*size);
    for (int strs = 0; strs < *size; ++strs) {
        fscanf(file,"%lg",&x[strs]);
    }
    return x;
}
