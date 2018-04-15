#ifndef _Util_H
#define _Util_H

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_blas.h"
#include "Qubit.h"

#define EPSILON 1E-5
#define BOOL_AS_STRING(b) ((b)? "True" : "False")

int print_with_indent(int indent, char* string, ...);
int print_double(double value);
int print_complex(gsl_complex value);
int print_matrix_complex(gsl_matrix_complex* m);
int print_vector_complex(gsl_vector_complex* v);
int print_qubit(Qubit* qubit);
int print_qupair(Qupair* qupair);
int print_qureg(Qureg* qureg);
gsl_vector_complex* matrix_vector_complex_mul(gsl_matrix_complex* m, gsl_vector_complex* v);
gsl_matrix_complex* matrix_matrix_complex_mul(gsl_matrix_complex* m1, gsl_matrix_complex* m2);
bool complex_is_zero(gsl_complex v);
bool complex_equal(gsl_complex a, gsl_complex b);
bool double_is_zero(double v);
bool double_equal(double a, double b);
int vector_complex_positions_swap(gsl_vector_complex* v, int sigs_num, int sig1, int sig2);

#endif
