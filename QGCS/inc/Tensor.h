#ifndef _Tensor_H
#define _Tensor_H

#include "Complex.h"

typedef struct _matrix {
    Complex* block;
    int size1;
    int size2;
} Matrix;

Matrix matrix_calloc(int size1, int size2);
Matrix matrix_from_array(Complex* mem, int size1, int size2);
Matrix matrix_from_matrix(Matrix m);
void matrix_free(Matrix m);
Complex matrix_get(Matrix m, int r, int c);
void matrix_set(Matrix m, int r, int c, Complex v);
void matrix_set_identity(Matrix a);
void matrix_swap(Matrix a, int r1, int c1, int r2, int c2);
Matrix matrix_add(Matrix a, Matrix b);
Matrix matrix_sub(Matrix a, Matrix b);
Matrix matrix_mul(Matrix a, Matrix b);
Matrix matrix_mul_scalar(Matrix a, Complex b);
Matrix matrix_trans(Matrix a);
Matrix matrix_conj(Matrix a);
Matrix matrix_dagger(Matrix a);

typedef struct _ket {
    Complex* block;
    int size;
} Ket;

typedef struct _bra {
    Complex* block;
    int size;
} Bra;

Ket ket_calloc(int size);
Bra bra_calloc(int size);
Ket ket_from_array(Complex* mem, int size);
Bra bra_from_array(Complex* mem, int size);
Ket ket_from_ket(Ket k);
Bra bra_from_bra(Bra b);
void ket_free(Ket k);
void bra_free(Bra b);
Bra ket_dagger(Ket a);
Ket bra_dagger(Bra a);
Complex ket_get(Ket a, int i);
void ket_set(Ket a, int i, Complex v);
Complex bra_get(Bra a, int i);
void bra_set(Bra a, int i, Complex v);
void ket_swap(Ket a, int i, int j);
void bra_swap(Bra a, int i, int j);
Ket matrix_mul_ket(Matrix a, Ket b);
Bra bra_mul_matrix(Bra a, Matrix b);
Complex inner_product(Bra a, Ket b);
Matrix outer_product(Ket a, Bra b);

Matrix matrix_Kronecker_product(Matrix a, Matrix b);
Ket ket_Kronecker_product(Ket a, Ket b);
int matrix_complex_swap(Matrix m, int r1, int c1, int r2, int c2);
int ket_complex_swap(Ket v, int x, int y);

#endif
