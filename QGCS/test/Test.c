#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include "Qubit.h"
#include "Util.h"
#include "Gate.h"
#include "Tensor.h"
#include "Algorithm.h"
#include "Feature.h"
#include "Memory.h"

int test_entanglement() {
    Qureg* qureg = allocate_qureg(3);
    for (int i = 0; i < 2; ++i) {
        apply_H(qureg->qubits[i]);
    }
    apply_CNOT(qureg->qubits, 3);
    print_qureg(qureg);

    apply_X(qureg->qubits[2]);
    print_qureg(qureg);

    apply_CNOT(qureg->qubits, 3);
    print_qureg(qureg);

    free_qureg(qureg);

    return 0;
}

int test_measurement() {
    Qureg* qureg = allocate_qureg(3);
    apply_H(qureg->qubits[0]);

    Qubit** temp = (Qubit**)pointer_memory_get(2);
    temp[0] = qureg->qubits[0];
    temp[1] = qureg->qubits[1];
    apply_CNOT(temp, 2);
    pointer_memory_return(temp, 2);

    print_qureg(qureg);

    apply_CNOT(qureg->qubits, 3);
    
    print_qureg(qureg);

    for (int i = 0; i < 3; ++i) {
        apply_X(qureg->qubits[i]);
    }
    
    print_qureg(qureg);

    apply_Z(qureg->qubits[1]);

    print_qureg(qureg);

    apply_PauliZ_M(qureg->qubits[1]);

    print_qureg(qureg);

    return 0;
}

int test_find_minimum() {
    int database_size = 128;
    int* database = int_memory_get(database_size);
    printf("QGCS\n\n");
    printf("Database size: %d\n", database_size);
    printf("Elements:\n");
    for (int i = 0; i < database_size; ++i) {
        database[i] = rand() % 100;
        printf("%d ", database[i]);
        if (!((i + 1) % 10)) {
            printf("\n");
        }
    }
    printf("\n");
    printf("\n");
    
    int min_index = find_minimum(database, database_size);
    
    int_memory_return(database, database_size);

    return 0;
}

int test_ket_positions_swap() {
    Ket v = ket_calloc(8);
    for (int i = 0; i < 8; ++i) {
        ket_set(v, i, complex_rect((double)i, 0));
    }
    print_ket(v);
    printf("\n");
    ket_positions_swap(v, 3, 0, 2);
    print_ket(v);
    printf("\n");

    return 0;
}

int test_possibility() {
    int repeat_times = 1000;
    int zero_count = 0;
    for (int i = 0; i < 1000; ++i) {
        Qureg* qureg = allocate_qureg(1);
        apply_H(qureg->qubits[0]);
        apply_PauliZ_M(qureg->qubits[0]);
        zero_count += qureg->qubits[0]->value == Zero ? 1 : 0;
        free_qureg(qureg);
    }
    printf("Repeat times = %d\n", repeat_times);
    printf("Zero count = %d\n", zero_count);

    return 0;
}

int main() {
    memory_init();
    unsigned int seed = (unsigned int)time(NULL); //1525627507 //1525627293
    srand(seed);
    test_find_minimum();
    //test_measurement();

    return 0;
}
