#include <stdio.h>
#include <stdarg.h>

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
    if(!qubit->entangled) {
        print_with_indent(indent, "state: (");
        print_complex(gsl_vector_complex_get(qubit->state, 0));
        printf(")|0> + (");
        print_complex(gsl_vector_complex_get(qubit->state, 1));
        printf(")|1>\n");
    } else {
        print_with_indent(indent, "<entangled>\n");
    }
    print_with_indent(indent, "measured: %s\n", BOOL_AS_STRING(qubit->measured));
    print_with_indent(indent, "value: %s\n", RESULT_AS_STRING(qubit->value));
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
        print_with_indent(indent, "|%*d>:    ", qupair->states_num / 10 + 1, i);
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

int print_complex(gsl_complex value) {
    printf("%.3f + %.3fi", GSL_REAL(value), GSL_IMAG(value));
    return 0;
}

int print_matrix_complex(gsl_matrix_complex* m) {
    for (size_t i = 0; i < m->size1; ++i) {
        printf("| ");
        for (size_t j = 0; j < m->size2; ++j) {
            print_complex(gsl_matrix_complex_get(m, i, j));
            printf("  ");
        }
        printf(" |\n");
    }
    return 0;
}

int print_vector_complex(gsl_vector_complex* v) {
    printf("| ");
    for (size_t i = 0; i < v->size; ++i) {
        print_complex(gsl_vector_complex_get(v, i));
        printf("  ");
    }
    printf(" |T");
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
    }
    printf("number of qupairs: %d\n", qureg->qupairs_num);
    printf("qupairs:\n");
    for(int i = 0; i < qureg->qupairs_num; ++i) {
        print_qupair_with_indent(qureg->qupairs[i], 1);
    }
    printf("\n");
    return 0;
}
