#ifndef _MATRIX_
#define _MATRIX_

#define NULL (void *)0;

typedef struct {
	double **values;
	unsigned int rows;
	unsigned int columns;
} matrix_t;


typedef struct {
	double *values;
	unsigned int length;
} vector_t;





#endif