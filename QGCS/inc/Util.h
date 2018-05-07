#ifndef _Util_H
#define _Util_H

#include "Qubit.h"

#define BOOL_AS_STRING(b) ((b) ? "True" : "False")

void print_with_indent(int indent, char* string, ...);
void print_double(double value);
void print_complex(Complex value);
void print_matrix(Matrix m);
void print_ket(Ket v);
void print_qubit(Qubit* qubit);
void print_qupair(Qupair* qupair);
void print_qureg(Qureg* qureg);

#endif
