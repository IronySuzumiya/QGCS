#include "Tensor.h"

gsl_matrix_complex* Kronecker_product_mm(gsl_matrix_complex* a, gsl_matrix_complex* b) {
    int size1 = a->size1 * b->size1;
    int size2 = a->size2 * b->size2;
    gsl_matrix_complex* result = gsl_matrix_complex_calloc(size1, size2);
    gsl_complex ca;
    gsl_complex cb;

    for (size_t i = 0; i < a->size1; ++i) {
        for (size_t j = 0; j < a->size2; ++j) {
            ca = gsl_matrix_complex_get(a, i, j);
            for (size_t m = 0; m < b->size1; ++m) {
                for (size_t n = 0; n < b->size2; ++n) {
                    cb = gsl_matrix_complex_get(b, m, n);
                    gsl_matrix_complex_set(result, i * b->size1 + m, j * b->size2 + n, gsl_complex_mul(ca, cb));
                }
            }
        }
    }
    return result;
}

gsl_vector_complex* Kronecker_product_vv(gsl_vector_complex* a, gsl_vector_complex* b) {
    int size = a->size * b->size;
    gsl_vector_complex* result = gsl_vector_complex_calloc(size);
    gsl_complex ca;
    gsl_complex cb;

    for (size_t i = 0; i < a->size; ++i) {
        ca = gsl_vector_complex_get(a, i);
        for (size_t j = 0; j < b->size; ++j) {
            cb = gsl_vector_complex_get(b, j);
            gsl_vector_complex_set(result, i * b->size + j, gsl_complex_mul(ca, cb));
        }
    }
    return result;
}
