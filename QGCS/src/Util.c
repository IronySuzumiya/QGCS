#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>

#include "Util.h"

#if _MSC_VER>=1900
#include "stdio.h" 
FILE* __cdecl __iob_func(unsigned i) {
    return __acrt_iob_func(i);
}
#endif /* _MSC_VER>=1900 */

static int print_qubit_with_indent(Qubit* qubit, int indent) {
    print_with_indent(indent, "index: %d\n", qubit->index);
    print_with_indent(indent, "qupair index: %d\n", qubit->qupair->index);
    if(!qubit->entangled && !qubit->measured) {
        print_with_indent(indent, "state: (");
        print_complex(gsl_vector_complex_get(qubit->state, 0));
        printf(")|0> + (");
        print_complex(gsl_vector_complex_get(qubit->state, 1));
        printf(")|1>\n");
    }
    else if(qubit->entangled) {
        print_with_indent(indent, "<entangled>\n");
    }
    else if (qubit->measured) {
        print_with_indent(indent, "<measured>\n");
        print_with_indent(indent, "value: %s\n", RESULT_AS_STRING(qubit->value));
    }
    return 0;
}

static int print_qupair_with_indent(Qupair* qupair, int indent) {
    print_with_indent(indent, "index: %d\n", qupair->index);
    print_with_indent(indent, "number of entangled qubits: %d\n", qupair->qubits_num);
    print_with_indent(indent, "indices of entangled qubits: ");
    for(int i = 0; i < qupair->qubits_num; ++i) {
        printf("%d ", qupair->qubits_indices[i]);
    }
    printf("(high to low)\n");
    print_with_indent(indent, "probability amplitude:\n");
    ++indent;
    for(int i = 0; i < qupair->states_num; ++i) {
        print_with_indent(indent, "|", qupair->states_num / 10 + 1, i);
        for (int j = 0; j < qupair->qubits_num; ++j) {
            printf("%d", (i >> (qupair->qubits_num - 1 - j)) & 0x1);
        }
        printf(">:    ");
        print_complex(gsl_vector_complex_get(qupair->state, i));
        printf("\n");
    }
    return 0;
}

int print_with_indent(int indent, char* string, ...) {
    va_list ap;
    fprintf(stdout, "%*s", indent * 4, "");
    va_start(ap, string);
    vfprintf(stdout, string, ap);
    va_end(ap);
    return 0;
}

int print_double(double value) {
    if (double_is_zero(value - round(value))) {
        printf("%d", (int)round(value));
    }
    else {
        printf("%.3lf", value);
    }
    return 0;
}

int print_complex(gsl_complex value) {
    double real = GSL_REAL(value);
    double imag = GSL_IMAG(value);
    if (complex_is_zero(value)) {
        printf("0");
    }
    else if(double_is_zero(real)) {
        print_double(imag);
        printf("i");
    }
    else if (double_is_zero(imag)) {
        print_double(real);
    }
    else {
        print_double(real);
        printf(" + ");
        print_double(imag);
        printf("i");
    }
    return 0;
}

int print_matrix_complex(gsl_matrix_complex* m) {
    for (size_t i = 0; i < m->size1; ++i) {
        printf(" ");
        for (size_t j = 0; j < m->size2; ++j) {
            print_complex(gsl_matrix_complex_get(m, i, j));
            printf("  ");
        }
        printf(" \n");
    }
    return 0;
}

int print_vector_complex(gsl_vector_complex* v) {
    printf("| ");
    for (size_t i = 0; i < v->size; ++i) {
        print_complex(gsl_vector_complex_get(v, i));
        printf("  ");
    }
    printf(" |");
    return 0;
}

int print_qubit(Qubit* qubit) {
    return print_qubit_with_indent(qubit, 0);
}

int print_qupair(Qupair* qupair) {
    return print_qupair_with_indent(qupair, 0);
}

int print_qureg(Qureg* qureg) {
    printf("number of qubits: %d\n", qureg->qubits_num);
    printf("qubits:\n");
    for(int i = 0; i < qureg->qubits_num; ++i) {
        print_qubit_with_indent(qureg->qubits[i], 1);
        printf("\n");
    }
    printf("number of qupairs: %d\n", qureg->qupairs_num);
    printf("qupairs:\n");
    for(int i = 0; i < qureg->qupairs_num; ++i) {
        print_qupair_with_indent(qureg->qupairs[i], 1);
        printf("\n");
    }
    printf("\n");
    return 0;
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

bool complex_is_zero(gsl_complex v) {
    return GSL_REAL(v) > -EPSILON && GSL_REAL(v) < EPSILON && GSL_IMAG(v) > -EPSILON && GSL_IMAG(v) < EPSILON;
}

bool complex_equal(gsl_complex a, gsl_complex b) {
    return complex_is_zero(gsl_complex_sub(a, b));
}

bool double_is_zero(double v) {
    return v > -EPSILON && v < EPSILON;
}

bool double_equal(double a, double b) {
    return double_is_zero(a - b);
}

int vector_complex_positions_swap(gsl_vector_complex* v, int sigs_num, int sig1, int sig2) {
    assert(sig1 < sigs_num && sig2 < sigs_num);
    assert(pow(2, sigs_num) == v->size);
    if (sig1 == sig2) {
        return 0;
    }
    bool* swapped = (bool*)malloc(sizeof(bool) * v->size);
    memset(swapped, 0, sizeof(bool) * v->size);
    for (int i = 0; i < (int)v->size; ++i) {
        // sigs are high-to-low
        if (swapped[i]) {
            continue;
        }
        int sig1_pos = (i >> (sigs_num - 1 - sig1)) & 0x1;
        int sig2_pos = (i >> (sigs_num - 1 - sig2)) & 0x1;
        if (sig1_pos == sig2_pos) {
            // don't need to swap
            swapped[i] = true;
            continue;
        }
        int target_i = 0;
        for (int j = 0; j < sigs_num; ++j) {
            if (j == sig1) {
                target_i |= sig2_pos << (sigs_num - 1 - j);
            }
            else if (j == sig2) {
                target_i |= sig1_pos << (sigs_num - 1 - j);
            }
            else {
                target_i |= i & (0x1 << (sigs_num - 1 - j));
            }
        }
        gsl_complex temp = gsl_vector_complex_get(v, i);
        gsl_vector_complex_set(v, i, gsl_vector_complex_get(v, target_i));
        gsl_vector_complex_set(v, target_i, temp);
        swapped[i] = true;
        swapped[target_i] = true;
    }
    free(swapped);
    return 0;
}
