#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "Qubit.h"

Qubit* allocate_qubit() {
    Qubit* qubit = (Qubit*)malloc(sizeof(Qubit));
    memset(qubit, 0, sizeof(Qubit));
    qubit->state = gsl_vector_complex_calloc(2);
    gsl_vector_complex_set(qubit->state, 0, gsl_complex_rect(1.0, 0));
    return qubit;
}

Qureg* allocate_qureg(int qubits_num) {
    Qureg* qureg = (Qureg*)malloc(sizeof(Qureg));
    qureg->qubits_num = qubits_num;
    qureg->qupairs_num = qubits_num;

    Qubit** qubits = (Qubit**)malloc(sizeof(Qubit*) * qureg->qubits_num);
    for(int i = 0; i < qureg->qubits_num; ++i) {
        qubits[i] = allocate_qubit();
        qubits[i]->index = i;
        qubits[i]->qureg = qureg;
    }
    qureg->qubits = qubits;

    Qupair** qupairs = (Qupair**)malloc(sizeof(Qupair*) * qureg->qupairs_num);
    for(int i = 0; i < qureg->qupairs_num; ++i) {
        qupairs[i] = (Qupair*)malloc(sizeof(Qupair));
        qupairs[i]->index = i;
        qupairs[i]->qureg = qureg;
        qupairs[i]->qubits_num = 1;
        qupairs[i]->qubits_indices = (int*)malloc(sizeof(int));
        qupairs[i]->qubits_indices[0] = i;
        qupairs[i]->states_num = 2;
        qupairs[i]->state = gsl_vector_complex_calloc(2);
        gsl_vector_complex_set(qupairs[i]->state, 0, gsl_complex_rect(1.0, 0));
        qubits[i]->qupair = qupairs[i];
    }
    qureg->qupairs = qupairs;

    return qureg;
}

int free_qubit(Qubit* qubit) {
    gsl_vector_complex_free(qubit->state);
    free(qubit);
    qubit = NULL;
    return 0;
}

int free_qupair(Qupair* qupair) {
    free(qupair->qubits_indices);
    gsl_vector_complex_free(qupair->state);
    free(qupair);
    return 0;
}

int free_qureg(Qureg* qureg) {
    for (int i = 0; i < qureg->qubits_num; ++i) {
        free_qubit(qureg->qubits[i]);
    }
    free(qureg->qubits);
    for (int i = 0; i < qureg->qupairs_num; ++i) {
        free_qupair(qureg->qupairs[i]);
    }
    free(qureg->qupairs);
    free(qureg);
    return 0;
}

int qubit_index_in_qupair(Qubit* qubit, Qupair* qupair) {
    assert(qubit->qureg == qupair->qureg);
    int index = -1;
    for (int i = 0; i < qupair->qubits_num; ++i) {
        if (qubit->index == qupair->qubits_indices[i]) {
            index = i;
            break;
        }
    }
    return index;
}
