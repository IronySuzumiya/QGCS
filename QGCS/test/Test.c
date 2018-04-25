#include <stdio.h>
#include <time.h>

#include "Qubit.h"
#include "Util.h"
#include "Gate.h"
#include "Tensor.h"
#include "Algorithm.h"

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

    apply_to_each(X.apply, qureg->qubits, 3);
    
    print_qureg(qureg);

    Z.apply(qureg->qubits[1]);

    print_qureg(qureg);

    PauliZ_M.apply(qureg->qubits[1]);

    print_qureg(qureg);

    return 0;
}

int test_find_minimum() {
    int database_size = 64;
    int* database = (int*)malloc(sizeof(int) * database_size);
    printf("Database is: ");
    for (int i = 0; i < database_size; ++i) {
        database[i] = rand() % 100;
        printf("%d ", database[i]);
    }
    printf("\n");

    int min_index = find_minimum(database, database_size);
    printf("Minimum's index is: %d\n", min_index);
    printf("Minimum's value is: %d\n", database[min_index]);

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

int test_possibility() {
    int repeat_times = 1000;
    int zero_count = 0;
    for (int i = 0; i < 1000; ++i) {
        Qureg* qureg = allocate_qureg(1);
        H.apply(qureg->qubits[0]);
        PauliZ_M.apply(qureg->qubits[0]);
        zero_count += qureg->qubits[0]->value == Zero ? 1 : 0;
        free_qureg(qureg);
    }
    printf("Repeat times = %d\n", repeat_times);
    printf("Zero count = %d\n", zero_count);

    return 0;
}

int test_H() {
    Qureg* qureg = allocate_qureg(3);
    apply_to_each(X.apply, qureg->qubits, qureg->qubits_num);
    CZ.apply(qureg->qubits, qureg->qubits_num);
    print_qureg(qureg);
    apply_to_each(X.apply, qureg->qubits, qureg->qubits_num);
    print_qureg(qureg);
    Qubit** control = (Qubit**)malloc(sizeof(Qubit*) * 2);
    for (int i = 0; i < qureg->qubits_num - 1; ++i) {
        control[i] = qureg->qubits[i];
    }
    apply_to_each(H.apply, control, qureg->qubits_num - 1);
    print_qureg(qureg);
    CNOT.apply(qureg->qubits, qureg->qubits_num);
    print_qureg(qureg);
    apply_to_each(H.apply, qureg->qubits, qureg->qubits_num);
    print_qureg(qureg);
    free_qureg(qureg);
    free(control);

    return 0;
}

int main() {
    gate_init(2333);
    //gate_init((unsigned int)time(NULL));
    //test_find_minimum();
    test_measurement();

    return 0;
}
