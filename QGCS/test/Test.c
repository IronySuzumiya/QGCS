#include "Qubit.h"
#include "Util.h"
#include "Gate.h"
#include "Tensor.h"
#include <stdio.h>

int main() {
    gate_init();

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
