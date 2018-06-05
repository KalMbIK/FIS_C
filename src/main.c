#include "suggested_work2.h"

void testCG(int* argc, char** argv[]);
void testGMRES(int* argc, char** argv[]);
void testMatvec(int* argc, char** argv[]);
double* readVec_from_vecfile(char* vecfile, int* size);

int main(int argc, char* argv[]) {
    testMatvec(&argc, &argv);
//    testGMRES(&argc, &argv);
    return 0;
}

void testCG(int* argc, char** argv[]){
    StorageFormat coo, csr;
    double *V = mm_read_mtx_to_COO((*argv)[3], &coo);
    formatConvertCOOtoCSR(&coo, &csr);
    size_t n = (size_t)csr.nx;

    double *x0 = newDoubleArray(n);
    double *b = newDoubleArray(n);
    double *x_star = ones(n);
    csr.matvec(&csr, V, x_star, b);

//    double *x = GMRES_full(&csr, V, x0, b, atoi((*argv)[1]));
    double  *x = CG_body(&csr, V, x0, b, atof((*argv)[2]));
    subtr(x,x_star,x_star,n);
    printDoubleArray(x, n);
    printf("||error||=%lg\n",getNorm(x_star,n));

    formatFree(&coo);
    formatFree(&csr);
    free(V);
    free(x0);
    free(x);
    free(x_star);
    free(b);
}

void testGMRES(int* argc, char** argv[]){
    StorageFormat coo, csr;
    double *V = mm_read_mtx_to_COO((*argv)[3], &coo);
    formatConvertCOOtoCSR(&coo, &csr);
    size_t n = (size_t)csr.nx;

    double *x0 = newDoubleArray(n);
    double *b = newDoubleArray(n);
    double *x_star = ones(n);
    csr.matvec(&csr, V, x_star, b);

//    double *x = GMRES_full(&csr, V, x0, b, atoi((*argv)[1]));
    double  *x = GMRES_restarted(&csr, V, x0, b, atoi((*argv)[1]), atof((*argv)[2]));
    subtr(x,x_star,x_star,n);
    printDoubleArray(x, n);
    printf("||error||=%lg\n",getNorm(x_star,n));

    formatFree(&coo);
    formatFree(&csr);
    free(V);
    free(x0);
    free(x);
    free(x_star);
    free(b);
}

void vectorOfIntsToFile(FILE *filePointer, int *vector, int size){
    fprintf(filePointer,"%d",vector[0]);
    for (int i = 1; i < size; i++){
        fprintf(filePointer,",%d",vector[i]);
    }
    fprintf(filePointer,"\n");
}

void testMatvec(int* argc, char** argv[]){
    StorageFormat coo, csr;
    double *V = mm_read_mtx_to_COO((*argv)[3], &coo);
    int vecSize;
    double *b2 = readVec_from_vecfile((*argv)[4],&vecSize);
    formatConvertCOOtoCSR(&coo, &csr);

    int n = csr.nx;
    FILE *kotu = fopen("kotocheck.csv","w");
    vectorOfIntsToFile(kotu,csr.I,n+1); // kotu
    fclose(kotu);

    double *b1 = newDoubleArray(n);
    double *units = newDoubleArray(n);

    for (int i = 0; i < n; ++i) {
        units[i]=i;
    }

    csr.matvec(&csr, V, units, b1);
//    coo.matvec(&coo, V, units, b1);

    printDoubleArray(b1,n);
    subtr(b1,b2,b1,n);
    printf("Norm(err) = %lg\n",getNorm(b1,n));

    formatFree(&coo);
    formatFree(&csr);
    free(V);
    free(b1);
    free(units);
    free(b2);
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
