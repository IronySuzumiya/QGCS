#ifndef _Gate_H
#define _Gate_H

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_blas.h"
#include "Qubit.h"

typedef void(*gate_apply)(Qubit*);
typedef void(*gate_with_1_parameter_apply)(Qubit*, double);
typedef void(*controlled_gate_apply)(Qubit**, int);
typedef void(*controlled_on_int_gate_apply)(Qubit**, int, int);
typedef void(*controlled_gate_with_1_parameter_apply)(Qubit**, int, double);
typedef void(*controlled_on_int_gate_with_1_parameter_apply)(Qubit**, int, int, double);
typedef void(*measurement_apply)(Qubit*);

typedef struct _gate {
    gate_apply apply;
    gate_apply apply_dagger;
    gsl_matrix_complex* matrix;
    gsl_matrix_complex* dagger_matrix;
} Gate;

typedef struct _gate_with_1_parameter {
    gate_with_1_parameter_apply apply;
    gate_with_1_parameter_apply apply_dagger;
} Gate1;

typedef struct _controlled_gate {
    controlled_gate_apply apply;
    controlled_gate_apply apply_dagger;
    controlled_on_int_gate_apply apply_on_int;
    controlled_on_int_gate_apply apply_on_int_dagger;
} ControlledGate;

typedef struct _controlled_gate_with_1_parameter {
    controlled_gate_with_1_parameter_apply apply;
    controlled_gate_with_1_parameter_apply apply_dagger;
    controlled_on_int_gate_with_1_parameter_apply apply_on_int;
    controlled_on_int_gate_with_1_parameter_apply apply_on_int_dagger;
} ControlledGate1;

typedef struct _measurement {
    measurement_apply apply;
    gsl_matrix_complex** operators;
} Measurement;

extern Gate H;
extern Gate X;
extern Gate1 R;
extern Gate Z;
extern ControlledGate CNOT;
extern ControlledGate1 CR;
extern ControlledGate CZ;
extern Measurement PauliZ_M;

int check_entanglement(Qupair* qupair);
int apply_to_each(gate_apply apply, Qubit** qubits, int qubits_num);
int apply_to_each_reverse(gate_apply apply, Qubit** qubits, int qubits_num);
int gate_init(unsigned int seed);

#endif
