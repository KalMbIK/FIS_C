//
// Created by pusheen on 04.01.2018.
//

#include "myMatrices.h"


TYPE** newMatrix(size_t hSize, size_t vSize){
    TYPE** res = (TYPE**)malloc(vSize * sizeof(TYPE*));
    TYPE* t = (TYPE*)calloc(hSize*vSize, sizeof(TYPE));
    for (int i = 0, j = 0; i < hSize*vSize; i+=hSize, j++)
        res[j]=t+i;
    return res;
}

void freeMatrix(TYPE **a){
    double *arr = *a;
    free(arr);
    free(a);
}

void printMatrix(TYPE **a, size_t hSize, size_t vSize){
    printf("[");
    for (int i = 0; i < vSize; i++)
        printArray(a[i],hSize);
    printf("]\n");
}

TYPE** newMatrixTransposed(TYPE **a, size_t hSize, size_t vSize){
    double **ta = newMatrix(vSize,hSize);
    transposeMatrix(a,ta,hSize,vSize);
    return ta;
}

void transposeMatrix(TYPE **a, TYPE **ta, size_t hSize, size_t vSize){
    for (int i = 0; i < hSize; i++){
        for (int j = 0; j < vSize; j++){
            ta[i][j]=a[j][i];
        }
    }
}

void diagMatrixVector(TYPE *a, TYPE *b, TYPE *res, size_t size){
    for(int i = 0; i < size; i++)
        res[i]=a[i]*b[i];
}
void matrixVector(TYPE **a, TYPE *b, TYPE *res, size_t hSize, size_t vSize) {
    for (int i = 0; i < vSize; i++){
        res[i] = dot(a[i],b,hSize);
    }
}

void matrixReduce(TYPE (*vFun)(TYPE *, size_t), TYPE **vectors, size_t size, TYPE *res, int numElements){
    for (int i = 0; i < numElements; i++)
        res[i] = vFun(vectors[i],size);
}