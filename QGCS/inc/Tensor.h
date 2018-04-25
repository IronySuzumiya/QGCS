#ifndef _Tensor_H
#define _Tensor_H

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_blas.h"

gsl_matrix_complex* Kronecker_product_mm(gsl_matrix_complex* a, gsl_matrix_complex* b);
gsl_vector_complex* Kronecker_product_vv(gsl_vector_complex* a, gsl_vector_complex* b);
gsl_vector_complex* matrix_vector_complex_mul(gsl_matrix_complex* m, gsl_vector_complex* v);
gsl_matrix_complex* matrix_matrix_complex_mul(gsl_matrix_complex* m1, gsl_matrix_complex* m2);
int matrix_complex_swap(gsl_matrix_complex* m, int r1, int c1, int r2, int c2);
int vector_complex_swap(gsl_vector_complex* v, int x, int y);

#endif
