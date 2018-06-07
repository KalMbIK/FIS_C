//
// Created by Hell on 17.05.2018.
//

#include "suggested_work1.h"

//CONSTRUCTORS and DESTRUCTOR
void formatSetParams(StorageFormat *format, int nx, int ny, int nz){
    format->nz = nz;
    format->ny = ny;
    format->nx = nx;
//    formatPrint(format);
}

void formatNewCOO(StorageFormat *coo, int nx, int ny, int nz) {
    formatSetParams(coo, nx, ny, nz);
    coo->I = (int*)calloc(nz, sizeof(int));
    coo->J = (int*)calloc(nz, sizeof(int));
    coo->matvec = &matvecCOO;
}

void formatNewCSR(StorageFormat *csr, int nx, int ny, int nz) {
    formatSetParams(csr, nx, ny, nz);
    csr->I = (int*)calloc(ny+1, sizeof(int));
    csr->J = (int*)calloc(nz, sizeof(int));
    csr->matvec = &matvecCSR;
}

void formatNewCSC(StorageFormat *csc, int nx, int ny, int nz) {
    formatSetParams(csc, nx, ny, nz);
    csc->I = (int*)calloc(nx+1, sizeof(int));
    csc->J = (int*)calloc(nz, sizeof(int));
    csc->matvec = &matvecCSC;
}

//CONVERSION
void formatConvertCOOtoCSR(StorageFormat *coo, StorageFormat *csr){
    formatNewCSR(csr,coo->nx,coo->ny,coo->nz);
    memcpy(csr->J,coo->J,coo->nz*sizeof(int));
    csr->I[0]=0;
    int counter = 1, row = 0;
    for (int i = 1, j = 0; i < csr->nz; ++i) {
        if (coo->I[i]==row)
            counter++;
        else{
            csr->I[j+1]=csr->I[j]+counter;
            counter=1;
            row=coo->I[i];
            j++;
        }
    }
    csr->I[csr->ny]=csr->I[csr->ny-1]+counter;
    if (coo->matvec==&matvecCOO_symm)
        csr->matvec=&matvecCSR_symm;
    else
        csr->matvec=&matvecCSR;
}

void formatFree(StorageFormat *format){
    free(format->I);
    free(format->J);
    format->nx = format->ny = format->nz = 0;
}

//ALGORITHMS
//MatrixVector multiplications
void matvecCOO(StorageFormat *coo, double *V, double *x, double *y){
    for (int i = 0; i < coo->ny; ++i) {
        y[i] = 0;
    }
    for (int i = 0; i < coo->nz; ++i) {
        y[coo->I[i]]+=V[i]*x[coo->J[i]];
    }
}

void matvecCSR(StorageFormat *csr, double *V, double *x, double *y) {
    for (int i = 0; i < csr->ny; ++i) {
        y[i] = 0;
        int i1 = csr->I[i];
        int i2 = csr->I[i+1];
        for (int j = i1; j < i2; ++j) {
            y[i] += V[j]*x[csr->J[j]];
        }
    }
}

void matvecCSC(StorageFormat *csc, double *V, double *x, double *y) {
    for (int i = 0; i < csc->nx; ++i) {
        int i1 = csc->I[i];
        int i2 = csc->I[i+1];
        for (int j = i1; j < i2; ++j) {
            y[csc->J[j]] += V[j]*x[i];
        }
    }
}

void matvecCOO_symm(StorageFormat *coo, double *V, double *x, double *y){
    for (int i = 0; i < coo->nx; ++i) {
        y[i] = 0;
    }
    for (int i = 0; i < coo->nz; ++i) {
        y[coo->I[i]]+=V[i]*x[coo->J[i]];
        if (coo->I[i] != coo->J[i])
            y[coo->J[i]]+=V[i]*x[coo->I[i]];
    }
}

void matvecCSR_symm(StorageFormat *csr, double *V, double *x, double *y){
    for (int i = 0; i < csr->ny; ++i) {
        y[i] = 0;
    }
    for (int i = 0; i < csr->ny; ++i) {
        int i1 = csr->I[i];
        int i2 = csr->I[i+1];
        for (int j = i1; j < i2; ++j) {
            y[i] += V[j]*x[csr->J[j]];
            if (i != csr->J[j])
                y[csr->J[j]] += V[j]*x[i];
        }
    }
}


void matvecCSR_leftSeidel(StorageFormat *csr, double *V, double *x, double *y){

}

//Matrix norms
double normFrobenius(StorageFormat *fmt, double*V){
    double norm = 0;
    for (int i = 0; i < fmt->nz; ++i) {
        norm += V[i]*V[i];
    }
    return sqrt(norm);
}

//OTHER USEFUL METHODS
//Comparator for sorting the array of MatrixElements
int cmpElements(const void * in1, const void * in2){
    MatrixElement *a = (MatrixElement*)in1;
    MatrixElement *b = (MatrixElement*)in2;
    int i = a->i - b->i, j = a->j - b->j;
    if (i < 0)
        return -1;
    if (i == 0){
        if (j < 0)
            return -1;
        else return 1;
    }
    return 1;
}

double* mm_read_mtx_to_COO(char *path, StorageFormat *coo) {
    int nx, ny, nz;
    MM_typecode matcode;
    FILE *file = fopen(path, "r");

    //Read the file into the array of structures and sort the resulting array
    mm_read_banner(file, &matcode);
    mm_read_mtx_crd_size(file, &ny, &nx, &nz);
    MatrixElement * mtrx = (MatrixElement*)calloc(nz, sizeof(MatrixElement));
    for (int i = 0; i < nz; i++) {
        fscanf(file, "%d %d %lg\n", &mtrx[i].i, &mtrx[i].j, &mtrx[i].v);
        mtrx[i].i--;
        mtrx[i].j--;
    }
    qsort(mtrx, nz, sizeof(MatrixElement), &cmpElements);

    //Save the data in the right format
    formatNewCOO(coo, nx, ny, nz);
    if (mm_is_symmetric(matcode)){
        coo->matvec = &matvecCOO_symm;
    } else
        coo->matvec = &matvecCOO;
    double *V = newDoubleArray(nz);
    for (int i = 0; i < nz; i++) {
        coo->I[i] = mtrx[i].i;
        coo->J[i] = mtrx[i].j;
        V[i] = mtrx[i].v;
    }

    free(mtrx);
    fclose(file);
    return V;
}

int* newIntArray(size_t size){
    return (int*)calloc(size,sizeof(double));
}

double* newDoubleArray(int size){
    return (double*)calloc(size,sizeof(double));
}

//PRINT
void formatPrint(StorageFormat *format){
    printf("Ny=%d, Nx=%d, Nz(nmbr of nonzero entries)=%d\n", format->ny, format->nx, format->nz);
}

void printIntArray(int *a, size_t size){
    for (int i = 0; i < size; i++){
        printf("%d ", a[i]);
    }
    printf("\n");
}

void printDoubleArray(double *a, size_t size){
    for (int i = 0; i < size; i++){
        printf("%lg ", a[i]);
    }
    printf("\n");
}