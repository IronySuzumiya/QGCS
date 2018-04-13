#ifndef _Qubit_H
#define _Qubit_H

#include <stdbool.h>

#include "gsl/gsl_complex_math.h"
#include "gsl/gsl_vector.h"

typedef enum _result {
    Unknown,
    Zero,
    One
} Result;

#define RESULT_AS_STRING(result) ((result) == Unknown? "Unknown" : (result) == Zero? "Zero" : "One")
#define RESULT_AS_INT(result) (assert((result) != Unknown), (result) == Zero ? 0 : 1)

typedef struct _qubit {
    int index;
    struct _qureg* qureg;
    struct _qupair* qupair;
    gsl_vector_complex* state;
    bool entangled;
    bool measured;
    Result value;
} Qubit;

typedef struct _qupair {
    int index;
    struct _qureg* qureg;
    int qubits_num;
    int* qubits_indices; // high to low
    int states_num;
    gsl_vector_complex* state;
} Qupair;

typedef struct _qureg {
    int qubits_num;
    int qupairs_num;
    Qubit** qubits;
    Qupair** qupairs;
} Qureg;

Qubit* allocate_qubit(void);
Qureg* allocate_qureg(int qubits_num);
// not complete delete
int free_qubit(Qubit* qubit);
// not complete delete
int free_qupair(Qupair* qupair);
int free_qureg(Qureg* qureg);
int qubit_index_in_qupair(Qubit* qubit, Qupair* qupair);

#endif
