#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <math.h>

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

    qureg->qubits = (Qubit**)malloc(sizeof(Qubit*) * qureg->qubits_num);
    for(int i = 0; i < qureg->qubits_num; ++i) {
        qureg->qubits[i] = allocate_qubit();
        qureg->qubits[i]->index = i;
        qureg->qubits[i]->qureg = qureg;
    }

    qureg->qupairs = (Qupair**)malloc(sizeof(Qupair*) * qureg->qupairs_num);
    for(int i = 0; i < qureg->qupairs_num; ++i) {
        qureg->qupairs[i] = (Qupair*)malloc(sizeof(Qupair));
        initialize_qupair_with_single_qubit_default(qureg->qupairs[i], qureg, i, qureg->qubits[i]);
    }

    return qureg;
}

int initialize_qupair_with_single_qubit(Qupair* qupair, Qureg* qureg, int index, Qubit* qubit, gsl_complex alpha, gsl_complex beta) {
    assert(qureg->qupairs[index] == qupair);
    assert(qubit->qureg == qureg);
    qupair->index = index;
    qupair->qureg = qureg;
    qupair->qubits_num = 1;
    qupair->qubits_indices = (int*)malloc(sizeof(int));
    qupair->qubits_indices[0] = qubit->index;
    qupair->states_num = 2;
    qupair->state = gsl_vector_complex_calloc(2);
    gsl_vector_complex_set(qupair->state, 0, alpha);
    gsl_vector_complex_set(qupair->state, 1, beta);
    qubit->qupair = qupair;

    return 0;
}

int initialize_qupair_with_single_qubit_default(Qupair* qupair, Qureg* qureg, int index, Qubit* qubit) {
    return initialize_qupair_with_single_qubit(qupair, qureg, index, qubit, gsl_complex_rect(1.0, 0), gsl_complex_rect(0, 0));
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

int qupair_set_probamp(Qupair* qupair, int index, gsl_complex value) {
    assert(index < qupair->states_num);
    gsl_vector_complex_set(qupair->state, index, value);
    return 0;
}

int qubit_set_probamp(Qubit* qubit, gsl_complex alpha, gsl_complex beta) {
    assert(!qubit->entangled);
    assert(!qubit->measured);
    gsl_vector_complex_set(qubit->state, 0, alpha);
    gsl_vector_complex_set(qubit->state, 1, beta);
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

int results_as_int(Qubit** qubits, int qubits_num) {
#if _QGCS_DEBUG
    for (int i = 0; i < qubits_num; ++i) {
        assert(qubits[i]->measured);
    }
#endif
    int result = 0;
    for (int i = 0; i < qubits_num; ++i) {
        if (qubits[i]->value == One) {
            result += (int)pow(2, qubits_num - 1 - i);
        }
    }
    return result;
}
