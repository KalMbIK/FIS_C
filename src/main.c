#include "suggested_work2.h"

int main(int* argc, char* argv[]) {
    char *path001 = "../data/fidap001.mtx";
    char *orsirr = "../data/orsirr_1.mtx";
    char *s3rmt3m3 = "../data/s3rmt3m3.mtx";
    char *myMatrixpath = "../data/myMatrix.mtx";
    char *mySymmMatrixpath = "../data/mySymmMatrix.mtx";
    StorageFormat coo, csr;
    double *V = mm_read_mtx_to_COO(orsirr, &coo);
    formatConvertCOOtoCSR(&coo, &csr);

    double *x0 = newDoubleArray(csr.nx);
    double *x_star = newDoubleArray(csr.nx);
    double *b = newDoubleArray(csr.nx);

    for (int i = 0; i < csr.nx; ++i) {
        x_star[i]=1;
    }

    csr.matvec(&csr, V, x_star, b);
    double  *x = restartedGMRES(&csr, V, b, x0, atoi(argv[1]), atof(argv[2]));
    subtr(x,x_star,x_star,csr.nx);
    printf("||error||=%lg\n",getNorm(x_star,csr.nx));
//    printDoubleArray(x, csr.nx);
//    csr.matvec(&csr, V, x, x_star);
//    subtr(b,x_star,x_star,csr.nx);
//    printf("%lg\n",getNorm(x_star,csr.nx));

    formatFree(&coo);
    formatFree(&csr);
    free(path001);
    free(orsirr);
    free(s3rmt3m3);
    free(myMatrixpath);
    free(mySymmMatrixpath);
    free(V);
    free(x0);
    free(x);
    free(x_star);
    free(b);
    return 0;
}