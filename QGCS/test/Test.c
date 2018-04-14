#include "Qubit.h"
#include "Util.h"
#include "Gate.h"
#include "Tensor.h"
#include <stdio.h>
#include <time.h>

int test_entanglement() {
    Qureg* qureg = allocate_qureg(3);
    for (int i = 0; i < 2; ++i) {
        H.apply(qureg->qubits[i]);
    }
    CNOT.apply(qureg->qubits, 3);
    print_qureg(qureg);

    X.apply(qureg->qubits[2]);
    print_qureg(qureg);

    CNOT.apply(qureg->qubits, 3);
    print_qureg(qureg);

    free_qureg(qureg);

    return 0;
}

int test_measurement() {
    Qureg* qureg = allocate_qureg(3);
    H.apply(qureg->qubits[0]);

    Qubit** temp = (Qubit**)malloc(sizeof(Qubit*) * 2);
    temp[0] = qureg->qubits[0];
    temp[1] = qureg->qubits[1];
    CNOT.apply(temp, 2);
    free(temp);

    print_qureg(qureg);

    CNOT.apply(qureg->qubits, 3);
    
    print_qureg(qureg);

    X.apply(qureg->qubits[2]);
    
    print_qureg(qureg);

    PauliZ_M.apply(qureg->qubits[1], 10000);

    print_qureg(qureg);

    return 0;
}

int test_vector_complex_positions_swap() {
    gsl_vector_complex* v = gsl_vector_complex_calloc(8);
    for (int i = 0; i < 8; ++i) {
        gsl_vector_complex_set(v, i, gsl_complex_rect(i, 0));
    }
    print_vector_complex(v);
    printf("\n");
    vector_complex_positions_swap(v, 3, 0, 2);
    print_vector_complex(v);
    printf("\n");

    return 0;
}

int main() {
    gate_init();

    test_measurement();

    return 0;
}
