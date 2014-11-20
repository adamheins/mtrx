#ifndef _MTRX_
#define _MTRX_



typedef struct {
	double **values;
	unsigned int rows;
	unsigned int columns;
} matrix_t;


typedef struct {
	double *values;
	unsigned int length;
} vector_t;


/**************************************************************************
 * Vectors                                                                *
 **************************************************************************/

vector_t *vctr_empty(unsigned int);
vector_t *vctr_zeros(unsigned int);
vector_t *vctr_ones(unsigned int);

void vctr_destroy(vector_t *);

vector_t *vctr_copy(vector_t *);
void vctr_print(vector_t *);

int vctr_eq_len(vector_t *, vector_t *);


double vctr_dot_prod(vector_t *, vector_t *);
vector_t *vctr_cross_prod(vector_t *, vector_t *);
double vctr_mag(vector_t *);


/**************************************************************************
* Matrices                                                                *
**************************************************************************/

// Initializers.
matrix_t *mtrx_zeros(unsigned int, unsigned int);
matrix_t *mtrx_ones(unsigned int, unsigned int);
matrix_t *mtrx_id(unsigned int);
matrix_t *mtrx_diag(vector_t *);

// Destructor.
void mtrx_destroy(matrix_t *);

matrix_t *mtrx_copy(matrix_t *);
void mtrx_print(matrix_t *);

int mtrx_is_sqr(matrix_t *);
int mtrx_eq_dim(matrix_t *, matrix_t *);

double mtrx_max(matrix_t *);
double mtrx_min(matrix_t *);


// Standard matrix arithmetic.
matrix_t *mtrx_mult(matrix_t *, matrix_t *);
matrix_t *mtrx_scale(matrix_t *, double);
matrix_t *mtrx_add(matrix_t *, matrix_t *);
matrix_t *mtrx_subtract(matrix_t *, matrix_t *);

vector_t *mtrx_mult_vctr(matrix_t *, vector_t *);

// Pointwise arithmetic.
matrix_t *mtrx_pw_mult(matrix_t *, matrix_t *);
matrix_t *mtrx_pw_div(matrix_t *, matrix_t *);
matrix_t *mtrx_pw_pow(matrix_t *, matrix_t *);

// Matrix manipulation.
void mtrx_row_swap(matrix_t *, unsigned int, unsigned int);
void mtrx_col_swap(matrix_t *, unsigned int, unsigned int);
void mtrx_scale_row(matrix_t *, unsigned int, double);
void mtrx_scale_col(matrix_t *, unsigned int, double);

// Row and column accessors and mutators.
vector_t *mtrx_get_row(matrix_t *, unsigned int);
vector_t *mtrx_get_col(matrix_t *, unsigned int);
void mtrx_set_col(matrix_t *, vector_t *, unsigned int);
void mtrx_set_row(matrix_t *, vector_t *, unsigned int);

matrix_t *mtrx_submatrix(matrix_t *, array_t *, array_t *);

matrix_t *mtrx_transpose(matrix_t *);
double mtrx_det(matrix_t *);
matrix_t *mtrx_inv(matrix_t *);

// Solves a linear system.
vector_t *mtrx_solve(matrix_t *, vector_t *);


#endif
