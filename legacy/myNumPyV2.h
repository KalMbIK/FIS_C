//
// Updated by pusheen on 17.12.17.
//

#ifndef NMPDE_C_MYGRID_H
#define NMPDE_C_MYGRID_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <stdarg.h>

#ifndef TYPE
#define TYPE double
#endif //TYPE

TYPE* newArray(size_t size);
void printArray(TYPE *a, size_t size);

//void arrayOfVectorsToCsv(char *filename, TYPE **vectors, int size, unsigned num);
//void vectorsToCsv(char *filename, int size, unsigned num, ...);

void vectorToFile(FILE *filePointer, TYPE *vector, size_t size);
void vectorsToCsv(char *filename, int size, unsigned num, ...);

double dot(TYPE *a, TYPE *b, size_t size);
double getNorm(TYPE *a, size_t size);
//TYPE reduce(TYPE (*f)(TYPE), TYPE *gr, int size);
TYPE getMaxElement(TYPE *a, size_t size);
TYPE getMinElement(TYPE *a, size_t size);

void add(TYPE *a, TYPE *b, TYPE *res, size_t size);
void subtr(TYPE *a, TYPE *b, TYPE *res, size_t size);

void constMult(TYPE *a, TYPE alpha, TYPE *res, size_t size);

void linspace(TYPE a, TYPE b, TYPE *gr, size_t size);
void vectorize(TYPE (*f)(TYPE), TYPE *gr, TYPE *res, size_t size);
void vectorizeOfArrayInPoint(TYPE (*f)(TYPE, TYPE), TYPE *gr, TYPE *res, size_t size, TYPE t);

void vectorize2D(TYPE (*f)(TYPE, TYPE), TYPE *a1, TYPE *a2, TYPE *res, size_t size);

#endif //NMPDE_C_MYGRID_H
