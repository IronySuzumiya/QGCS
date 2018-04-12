#ifndef _Util_H
#define _Util_H

#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "Qubit.h"

#define BOOL_AS_STRING(b) ((b)? "True" : "False")

int print_with_indent(int indent, char* string, ...);
int print_complex(gsl_complex value);
int print_matrix_complex(gsl_matrix_complex* m);
int print_vector_complex(gsl_vector_complex* v);
int print_qubit(Qubit* qubit);
int print_qupair(Qupair* qupair);
int print_qureg(Qureg* qureg);

#endif
