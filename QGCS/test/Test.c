#include "Qubit.h"
#include "Util.h"
#include "Gate.h"
#include "Tensor.h"
#include <stdio.h>

int main() {
    gate_init();

    Qureg* qureg = allocate_qureg(5);
    print_qureg(qureg);

    H.apply(qureg->qubits[0]);
    print_qubit(qureg->qubits[0]);
    print_qupair(qureg->qubits[0]->qupair);

    gsl_matrix_complex* identity = gsl_matrix_complex_calloc(2, 2);
    gsl_matrix_complex_set_identity(identity);
    gsl_matrix_complex* IH = Kronecker_product(identity, H.matrix);

    printf("Kronecker Product of Hadamard Matrix and Identity Matrix:\n");
    print_matrix_complex(IH);

    //RESULT_AS_INT(qureg->qubits[0]->value);

    return 0;
}
