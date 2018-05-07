#include <math.h>
#include <assert.h>

#include "Qubit.h"
#include "Memory.h"

Qubit* allocate_qubit() {
    Qubit* qubit = qubit_memory_get(1);
    qubit->index = 0;
    qubit->qureg = 0;
    qubit->qupair = 0;
    qubit->state = ket_calloc(2);
    ket_set(qubit->state, 0, complex_rect(1.0, 0));
    ket_set(qubit->state, 1, complex_rect(0, 0));
    qubit->entangled = 0;
    qubit->measured = 0;
    qubit->value = Unknown;

    return qubit;
}

Qureg* allocate_qureg(int qubits_num) {
    Qureg* qureg = qureg_memory_get(1);
    qureg->qubits_num = qubits_num;
    qureg->qupairs_num = qubits_num;

    qureg->qubits = (Qubit**)pointer_memory_get(qubits_num);
    for(int i = 0; i < qubits_num; ++i) {
        qureg->qubits[i] = allocate_qubit();
        qureg->qubits[i]->index = i;
        qureg->qubits[i]->qureg = qureg;
    }

    qureg->qupairs = (Qupair**)pointer_memory_get(qureg->qupairs_num);
    for(int i = 0; i < qureg->qupairs_num; ++i) {
        qureg->qupairs[i] = qupair_memory_get(1);
        initialize_qupair_with_single_qubit_default(qureg->qupairs[i], qureg, i, qureg->qubits[i]);
    }

    return qureg;
}

int initialize_qupair_with_single_qubit(Qupair* qupair, Qureg* qureg, int index, Qubit* qubit, Complex alpha, Complex beta) {
    assert(qureg->qupairs[index] == qupair);
    assert(qubit->qureg == qureg);
    qupair->index = index;
    qupair->qureg = qureg;
    qupair->qubits_num = 1;
    qupair->qubits_indices = int_memory_get(1);
    qupair->qubits_indices[0] = qubit->index;
    qupair->states_num = 2;
    qupair->state = ket_calloc(2);
    ket_set(qupair->state, 0, alpha);
    ket_set(qupair->state, 1, beta);
    qubit->qupair = qupair;

    return 0;
}

int initialize_qupair_with_single_qubit_default(Qupair* qupair, Qureg* qureg, int index, Qubit* qubit) {
    return initialize_qupair_with_single_qubit(qupair, qureg, index, qubit, complex_rect(1.0, 0), complex_rect(0, 0));
}

int free_qubit(Qubit* qubit) {
    ket_free(qubit->state);
    qubit_memory_return(qubit, 1);
    qubit = 0;
    return 0;
}

int free_qupair(Qupair* qupair) {
    int_memory_return(qupair->qubits_indices, qupair->qubits_num);
    ket_free(qupair->state);
    qupair_memory_return(qupair, 1);
    qupair = 0;
    return 0;
}

int free_qureg(Qureg* qureg) {
    for (int i = 0; i < qureg->qubits_num; ++i) {
        free_qubit(qureg->qubits[i]);
    }
    pointer_memory_return(qureg->qubits, qureg->qubits_num);
    qureg->qubits = 0;
    for (int i = 0; i < qureg->qupairs_num; ++i) {
        free_qupair(qureg->qupairs[i]);
    }
    pointer_memory_return(qureg->qupairs, qureg->qupairs_num);
    qureg->qupairs = 0;
    qureg_memory_return(qureg, 1);
    qureg = 0;
    return 0;
}

int qupair_set_probamp(Qupair* qupair, int index, Complex value) {
    assert(index < qupair->states_num);
    ket_set(qupair->state, index, value);
    return 0;
}

int qubit_set_probamp(Qubit* qubit, Complex alpha, Complex beta) {
    assert(!qubit->entangled);
    assert(!qubit->measured);
    ket_set(qubit->state, 0, alpha);
    ket_set(qubit->state, 1, beta);
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
    for (int i = 0; i < qubits_num; ++i) {
        assert(qubits[i]->measured);
    }
    int result = 0;
    for (int i = 0; i < qubits_num; ++i) {
        if (qubits[i]->value == One) {
            result += (int)pow(2, qubits_num - 1 - i);
        }
    }
    return result;
}

void check_qubit(Qubit* qubit) {
    double alpha2 = complex_norm(ket_get(qubit->state, 0));
    double beta2 = complex_norm(ket_get(qubit->state, 1));
    assert(double_equal(alpha2 + beta2, 1));
}

void check_qupair(Qupair* qupair) {
    double sum = 0;
    for (int j = 0; j < qupair->states_num; ++j) {
        sum += complex_norm(ket_get(qupair->state, j));
    }
    if (!double_equal(sum, 1)) {
        ;
    }
    //assert(double_equal(sum, 1));
}

void check_qureg(Qureg* qureg) {
    int qubits_num = qureg->qubits_num;
    Qubit** qubits = qureg->qubits;
    int* checked = int_memory_get(qubits_num);
    for (int i = 0; i < qubits_num; ++i) {
        checked[i] = 0;
    }
    for (int i = 0; i < qubits_num; ++i) {
        if (!checked[i]) {
            if (!qubits[i]->entangled) {
                check_qubit(qubits[i]);
                checked[i] = 1;
            }
            else {
                Qupair* qupair = qubits[i]->qupair;
                check_qupair(qupair);
                for (int j = 0; j < qupair->qubits_num; ++j) {
                    checked[qupair->qubits_indices[i]] = 1;
                }
            }
        }
    }
    int_memory_return(checked, qubits_num);
}
