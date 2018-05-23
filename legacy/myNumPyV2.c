//
// Updated by pusheen on 17.12.17.
//

#include "myNumPyV2.h"

TYPE* newArray(size_t size){
    return (TYPE*)calloc(size,sizeof(TYPE));
}

void printArray(TYPE *a, size_t size){
    for (int i = 0; i < size; i++){
        printf("%lg ", a[i]);
    }
    printf("\n");
}

//void arrayOfVectorsToCsv(char *filename, TYPE **vectors, int size, unsigned num){
//    FILE *fp;
//    fp=fopen(filename,"w+");
//    for (int i = 0; i < size; i++){
//        fprintf(fp,"%lg",vectors[0][i]);
//        for (int j = 1; j < num; j++){
//            fprintf(fp,",%lg",vectors[j][i]);
//        }
//        fprintf(fp,"\n");
//    }
//    fclose(fp);
//}
//
//void vectorsToCsv(char *filename, int size, unsigned num, ...){
//    va_list args;
//    va_start(args, num);
//    TYPE** vectors = (TYPE**)malloc(num*sizeof(TYPE*));
//    int ind = 0;
//    while (num--){
//        vectors[ind] = va_arg(args,TYPE*);
//        ind++;
//    }
//    va_end(args);
//    arrayOfVectorsToCsv(filename, vectors, size, ind);
//}

void vectorToFile(FILE *filePointer, TYPE *vector, size_t size){
    fprintf(filePointer,"%lg",vector[0]);
    for (int i = 1; i < size; i++){
        fprintf(filePointer,",%lg",vector[i]);
    }
    fprintf(filePointer,"\n");
}

void vectorsToCsv(char *filename, int size, unsigned num, ...){
    FILE *fp=fopen(filename,"w+");
    va_list args;
    va_start(args, num);
    while (num--){
        vectorToFile(fp, va_arg(args, TYPE *), size);
    }
    va_end(args);
    fclose(fp);
}

double dot(TYPE *a, TYPE *b, size_t size){
    double sum = 0;
    for (int i = 0; i < size; i++)
        sum += a[i]*b[i];
    return sum;
}

double getNorm(TYPE *a, size_t size){
    return sqrt(dot(a,a,size));
}

//TYPE reduce(TYPE (*f)(TYPE), TYPE *gr, int size){
//
//}

TYPE getMaxElement(TYPE *a, size_t size){
    double max = INT_MIN;
    for (int i = 0; i < size; i++){
        if (a[i] > max)
            max = a[i];
    }
    return max;
}

TYPE getMinElement(TYPE *a, size_t size){
    double min = INT_MAX;
    for (int i = 0; i < size; i++){
        if (a[i] < min)
            min = a[i];
    }
    return min;
}

void add(TYPE *a, TYPE *b, TYPE *res, size_t size){
    for (int i = 0; i < size; i++)
        res[i] = a[i]+b[i];
}

void subtr(TYPE *a, TYPE *b, TYPE *res, size_t size){
    for (int i = 0; i < size; i++)
        res[i] = a[i]-b[i];
}

void constMult(TYPE *a, TYPE alpha, TYPE *res, size_t size){
    for (int i = 0; i < size; i++)
        res[i]=alpha*a[i];
}

void linspace(TYPE a, TYPE b, TYPE *gr, size_t size) {
    TYPE step = (b-a)/(size-1);
    gr[0] = a;
    for (int i = 1; i < size; i++){
        gr[i]=gr[i-1]+step;
    }
}

void vectorize(TYPE (*f)(TYPE), TYPE *gr, TYPE *res, size_t size){
    for (int i = 0; i < size; i++){
        res[i]=f(gr[i]);
    }
}

void vectorizeOfArrayInPoint(TYPE (*f)(TYPE, TYPE), TYPE *gr, TYPE *res, size_t size, TYPE t){
    for (int i = 0; i < size; i++){
        res[i]=f(gr[i],t);
    }
}

void vectorize2D(TYPE (*f)(TYPE, TYPE), TYPE *a1, TYPE *a2, TYPE *res, size_t size){
    for (int i = 0; i < size; i++){
        res[i]=f(a1[i],a2[i]);
    }
}