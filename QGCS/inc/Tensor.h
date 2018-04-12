#ifndef _Tensor_H
#define _Tensor_H

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_complex_math.h"

gsl_matrix_complex* Kronecker_product_mm(gsl_matrix_complex* a, gsl_matrix_complex* b);
gsl_matrix_complex* Kronecker_product_vv(gsl_vector_complex* a, gsl_vector_complex* b);

#endif
