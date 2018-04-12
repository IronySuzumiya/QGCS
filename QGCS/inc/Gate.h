#ifndef _Gate_H
#define _Gate_H

#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_blas.h"
#include "Qubit.h"

typedef void(*gate_apply)(Qubit*);
typedef void(*controlled_gate_apply)(Qubit**, int);

typedef struct _gate {
    gate_apply apply;
    gate_apply apply_dagger;
    gsl_matrix_complex* matrix;
    gsl_matrix_complex* dagger_matrix;
} Gate;

typedef struct _controlled_gate {
    controlled_gate_apply apply;
    controlled_gate_apply apply_dagger;
} ControlledGate;

extern Gate H;
extern Gate X;
extern ControlledGate CNOT;

int gate_init(void);

#endif
