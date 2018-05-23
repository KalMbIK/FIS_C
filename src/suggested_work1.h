//
// Created by Hell on 17.05.2018.
//

#ifndef FIS_C_STORAGEFORMATSLIB_H
#define FIS_C_STORAGEFORMATSLIB_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "../lib/mmio.h"
#include "../legacy/myNumPyV2.h"

#define MATVEC_DEFAULT &matvecCSR

//With this data structure we could store the matrix as discussed in class
struct StorageFormat{
    void (*matvec)(struct StorageFormat*, double*, double*, double*);
    int nx, ny, nz;
    int *I, *J;
};
typedef struct StorageFormat StorageFormat;

//With that structure we can store the matrix as an array of such elements
struct MatrixElement{
    int i, j;
    double v;
};
typedef struct MatrixElement MatrixElement;

//CONSTRUCTORS and DESTRUCTOR
void formatSetParams(StorageFormat *format, int nx, int ny, int nz);
void formatNewCOO(StorageFormat *coo, int nx, int ny, int nz);
void formatNewCSR(StorageFormat *csr, int nx, int ny, int nz);
void formatNewCSC(StorageFormat *csc, int nx, int ny, int nz);
//CONVERSION
void formatConvertCOOtoCSR(StorageFormat *coo, StorageFormat *csr);
void formatFree(StorageFormat *format);

//ALGORITHMS
//MatrixVector multiplications
void matvecCSR(StorageFormat *csr, double *V, double *x, double *y);
void matvecCSC(StorageFormat *csc, double *V, double *x, double *y);
void matvecSymm(StorageFormat *fmt, double *V, double *x, double *y);

//Matrix norms
double normFrobenius(StorageFormat *fmt, double*V);

//OTHER USEFUL METHODS
//Comparator for sorting the array of MatrixElements
int cmpElements(const void * in1, const void * in2);
//Read the mtx file and store it as COO (Values array is returned)
double* mm_read_mtx_to_COO(char *path, StorageFormat *coo);

int* newIntArray(size_t size);
double* newDoubleArray(size_t size);
//PRINT
void formatPrint(StorageFormat *format);
void printIntArray(int *a, size_t size);
void printDoubleArray(double *a, size_t size);

#endif //FIS_C_STORAGEFORMATSLIB_H
