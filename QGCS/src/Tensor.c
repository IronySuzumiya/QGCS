#include <assert.h>

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

gsl_vector_complex* matrix_vector_complex_mul(gsl_matrix_complex* m, gsl_vector_complex* v) {
    assert(m->size2 == v->size);
    gsl_vector_complex* temp = gsl_vector_complex_calloc(v->size);
    gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(1.0, 0), m, v, gsl_complex_rect(0, 0), temp);

    return temp;
}

gsl_matrix_complex* matrix_matrix_complex_mul(gsl_matrix_complex* m1, gsl_matrix_complex* m2) {
    assert(m1->size2 == m2->size1);
    gsl_matrix_complex* temp = gsl_matrix_complex_calloc(m1->size1, m2->size2);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1.0, 0), m1, m2, gsl_complex_rect(0, 0), temp);

    return temp;
}

int matrix_complex_swap(gsl_matrix_complex* m, int r1, int c1, int r2, int c2) {
    gsl_complex temp = gsl_matrix_complex_get(m, r1, c1);
    gsl_matrix_complex_set(m, r1, c1, gsl_matrix_complex_get(m, r2, c2));
    gsl_matrix_complex_set(m, r2, c2, temp);

    return 0;
}

int vector_complex_swap(gsl_vector_complex* v, int x, int y) {
    gsl_complex temp = gsl_vector_complex_get(v, x);
    gsl_vector_complex_set(v, x, gsl_vector_complex_get(v, y));
    gsl_vector_complex_set(v, y, temp);

    return 0;
}
