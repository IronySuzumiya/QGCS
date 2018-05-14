#include <assert.h>

#include "Tensor.h"
#include "Memory.h"

Matrix matrix_calloc(int size1, int size2) {
    Complex* mem = complex_memory_get(size1 * size2);
    for (int i = 0; i < size1 * size2; ++i) {
        mem[i] = complex_rect(0, 0);
    }
    return (Matrix) { mem, size1, size2 };
}

Matrix matrix_from_array(Complex* mem, int size1, int size2) {
    return (Matrix) { mem, size1, size2 };
}

Matrix matrix_from_matrix(Matrix m) {
    Complex* mem = complex_memory_get(m.size1 * m.size2);
    for (int i = 0; i < m.size1; ++i) {
        for (int j = 0; j < m.size2; ++j) {
            mem[i * m.size2 + j] = matrix_get(m, i, j);
        }
    }
    return (Matrix) { mem, m.size1, m.size2 };
}

void matrix_free(Matrix m) {
    complex_memory_return(m.block, m.size1 * m.size2);
}

Complex matrix_get(Matrix m, int r, int c) {
    return m.block[r * m.size2 + c];
}

void matrix_set(Matrix m, int r, int c, Complex v) {
    m.block[r * m.size2 + c] = v;
}

void matrix_set_identity(Matrix a) {
    assert(a.size1 == a.size2);
    for (int i = 0; i < a.size1; ++i) {
        for (int j = 0; j < a.size2; ++j) {
            matrix_set(a, i, j, complex_rect(0, 0));
        }
        matrix_set(a, i, i, complex_rect(1, 0));
    }
}

void matrix_swap(Matrix a, int r1, int c1, int r2, int c2) {
    Complex temp = matrix_get(a, r1, c1);
    matrix_set(a, r1, c1, matrix_get(a, r2, c2));
    matrix_set(a, r2, c2, temp);
}

Matrix matrix_add(Matrix a, Matrix b) {
    assert(a.size1 == b.size1);
    assert(a.size2 == b.size2);
    Matrix c = matrix_calloc(a.size1, a.size2);
    for (int i = 0; i < c.size1; ++i) {
        for (int j = 0; j < c.size2; ++j) {
            matrix_set(c, i, j, complex_add(matrix_get(a, i, j), matrix_get(b, i, j)));
        }
    }
    return c;
}

Matrix matrix_sub(Matrix a, Matrix b) {
    assert(a.size1 == b.size1);
    assert(a.size2 == b.size2);
    Matrix c = matrix_calloc(a.size1, a.size2);
    for (int i = 0; i < c.size1; ++i) {
        for (int j = 0; j < c.size2; ++j) {
            matrix_set(c, i, j, complex_sub(matrix_get(a, i, j), matrix_get(b, i, j)));
        }
    }
    return c;
}

Matrix matrix_mul(Matrix a, Matrix b) {
    assert(a.size2 == b.size1);
    Matrix c = matrix_calloc(a.size1, b.size2);
    for (int i = 0; i < c.size1; ++i) {
        for (int j = 0; j < c.size2; ++j) {
            Complex temp = complex_rect(0, 0);
            for (int m = 0; m < a.size2; ++m) {
                temp = complex_add(temp, complex_mul(matrix_get(a, i, m), matrix_get(b, m, j)));
            }
            matrix_set(c, i, j, temp);
        }
    }
    return c;
}

Matrix matrix_mul_scalar(Matrix a, Complex b) {
    Matrix c = matrix_calloc(a.size1, a.size2);
    for (int i = 0; i < c.size1; ++i) {
        for (int j = 0; j < c.size2; ++j) {
            matrix_set(c, i, j, complex_mul(matrix_get(a, i, j), b));
        }
    }
    return c;
}

Matrix matrix_trans(Matrix a) {
    Matrix b = matrix_calloc(a.size2, a.size1);
    for (int i = 0; i < a.size1; ++i) {
        for (int j = i + 1; j < a.size2; ++j) {
            matrix_swap(b, i, j, j, i);
        }
    }
    return b;
}

Matrix matrix_conj(Matrix a) {
    Matrix b = matrix_from_matrix(a);
    for (int i = 0; i < b.size1; ++i) {
        for (int j = 0; j < b.size2; ++j) {
            matrix_set(b, i, j, complex_conj(matrix_get(b, i, j)));
        }
    }
    return b;
}

Matrix matrix_dagger(Matrix a) {
    Matrix b = matrix_calloc(a.size2, a.size1);
    Matrix c = matrix_trans(b);
    matrix_free(b);
    Matrix d = matrix_conj(c);
    matrix_free(c);
    return d;
}

Ket ket_calloc(int size) {
    Complex* mem = complex_memory_get(size);
    for (int i = 0; i < size; ++i) {
        mem[i] = complex_rect(0, 0);
    }
    return (Ket) { mem, size };
}

Bra bra_calloc(int size) {
    Complex* mem = complex_memory_get(size);
    for (int i = 0; i < size; ++i) {
        mem[i] = complex_rect(0, 0);
    }
    return (Bra) { mem, size };
}

Ket ket_from_array(Complex* mem, int size) {
    return (Ket) { mem, size };
}

Bra bra_from_array(Complex* mem, int size) {
    return (Bra) { mem, size };
}

Ket ket_from_ket(Ket k) {
    Complex* mem = complex_memory_get(k.size);
    for (int i = 0; i < k.size; ++i) {
        mem[i] = ket_get(k, i);
    }
    return (Ket) { mem, k.size };
}

Bra bra_from_bra(Bra b) {
    Complex* mem = complex_memory_get(b.size);
    for (int i = 0; i < b.size; ++i) {
        mem[i] = bra_get(b, i);
    }
    return (Bra) { mem, b.size };
}

void ket_free(Ket k) {
    complex_memory_return(k.block, k.size);
}

void bra_free(Bra b) {
    complex_memory_return(b.block, b.size);
}

Bra ket_dagger(Ket a) {
    Complex* mem = complex_memory_get(a.size);
    for (int i = 0; i < a.size; ++i) {
        mem[i] = complex_conj(ket_get(a, i));
    }
    return (Bra) { mem, a.size };
}

Ket bra_dagger(Bra a) {
    Complex* mem = complex_memory_get(a.size);
    for (int i = 0; i < a.size; ++i) {
        mem[i] = complex_conj(bra_get(a, i));
    }
    return (Ket) { mem, a.size };
}

Complex ket_get(Ket a, int i) {
    return a.block[i];
}

void ket_set(Ket a, int i, Complex v) {
    a.block[i] = v;
}

Complex bra_get(Bra a, int i) {
    return a.block[i];
}

void bra_set(Bra a, int i, Complex v) {
    a.block[i] = v;
}

void ket_swap(Ket a, int i, int j) {
    Complex temp = ket_get(a, i);
    ket_set(a, i, ket_get(a, j));
    ket_set(a, j, temp);
}

void bra_swap(Bra a, int i, int j) {
    Complex temp = bra_get(a, i);
    bra_set(a, i, bra_get(a, j));
    bra_set(a, j, temp);
}

Ket matrix_mul_ket(Matrix a, Ket b) {
    assert(a.size2 == b.size);
    Ket c = ket_calloc(a.size1);
    for (int i = 0; i < c.size; ++i) {
        Complex temp = complex_rect(0, 0);
        for (int j = 0; j < b.size; ++j) {
            temp = complex_add(temp, complex_mul(matrix_get(a, i, j), ket_get(b, j)));
        }
        ket_set(c, i, temp);
    }

    /*double sum = 0;
    for (int i = 0; i < c.size; ++i) {
        sum += complex_norm(ket_get(c, i));
    }
    assert(double_equal(sum, 1));*/

    return c;
}

Bra bra_mul_matrix(Bra a, Matrix b) {
    assert(a.size == b.size1);
    Bra c = bra_calloc(b.size2);
    for (int i = 0; i < c.size; ++i) {
        Complex temp = complex_rect(0, 0);
        for (int j = 0; j < a.size; ++j) {
            temp = complex_add(temp, complex_mul(bra_get(a, j), matrix_get(b, j, i)));
        }
        bra_set(c, i, temp);
    }
    return c;
}

Complex inner_product(Bra a, Ket b) {
    assert(a.size == b.size);
    Complex c = complex_rect(0, 0);
    for (int i = 0; i < a.size; ++i) {
        c = complex_add(c, complex_mul(bra_get(a, i), ket_get(b, i)));
    }
    return c;
}

Matrix outer_product(Ket a, Bra b) {
    Matrix c = matrix_calloc(a.size, b.size);
    for (int i = 0; i < a.size; ++i) {
        for (int j = 0; j < b.size; ++j) {
            matrix_set(c, i, j, complex_mul(ket_get(a, i), bra_get(b, j)));
        }
    }
    return c;
}

Matrix matrix_Kronecker_product(Matrix a, Matrix b) {
    int size1 = a.size1 * b.size1;
    int size2 = a.size2 * b.size2;
    Matrix result = matrix_calloc(size1, size2);
    Complex ca;
    Complex cb;

    for (int i = 0; i < a.size1; ++i) {
        for (int j = 0; j < a.size2; ++j) {
            ca = matrix_get(a, i, j);
            for (int m = 0; m < b.size1; ++m) {
                for (int n = 0; n < b.size2; ++n) {
                    cb = matrix_get(b, m, n);
                    matrix_set(result, i * b.size1 + m, j * b.size2 + n, complex_mul(ca, cb));
                }
            }
        }
    }

    return result;
}

Ket ket_Kronecker_product(Ket a, Ket b) {
    int size = a.size * b.size;
    Ket result = ket_calloc(size);
    Complex ca;
    Complex cb;

    for (int i = 0; i < a.size; ++i) {
        ca = ket_get(a, i);
        for (int j = 0; j < b.size; ++j) {
            cb = ket_get(b, j);
            ket_set(result, i * b.size + j, complex_mul(ca, cb));
        }
    }

    return result;
}
