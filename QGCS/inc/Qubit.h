#ifndef _Qubit_H
#define _Qubit_H

#include "Tensor.h"

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
    Ket state;
    int entangled;
    int measured;
    Result value;
} Qubit;

typedef struct _qupair {
    int index;
    struct _qureg* qureg;
    int qubits_num;
    int* qubits_indices; // high to low
    int states_num;
    Ket state;
} Qupair;

typedef struct _qureg {
    int qubits_num;
    int qupairs_num;
    Qubit** qubits;
    Qupair** qupairs;
} Qureg;

Qubit* allocate_qubit(void);
Qureg* allocate_qureg(int qubits_num);
int initialize_qupair_with_single_qubit(Qupair* qupair, Qureg* qureg, int index, Qubit* qubit, Complex alpha, Complex beta);
int initialize_qupair_with_single_qubit_default(Qupair* qupair, Qureg* qureg, int index, Qubit* qubit);
int free_qubit(Qubit* qubit);
int free_qupair(Qupair* qupair);
int free_qureg(Qureg* qureg);
int qupair_set_probamp(Qupair* qupair, int index, Complex value);
int qubit_set_probamp(Qubit* qubit, Complex alpha, Complex beta);
int qubit_index_in_qupair(Qubit* qubit, Qupair* qupair);
int results_as_int(Qubit** qubits, int qubits_num);
void check_qubit(Qubit* qubit);
void check_qupair(Qupair* qupair);
void check_qureg(Qureg* qureg);

#endif
