#include <math.h>
#include <assert.h>
#include <time.h>
#include <stdlib.h>

#include "Algorithm.h"
#include "Qubit.h"
#include "Gate.h"
#include "Util.h"
#include "Tensor.h"
#include "Feature.h"
#include "Const.h"
#include "Memory.h"

static double calculate_success_probability(int marked_count, int database_size, int iterations_num) {
    double a = asin(sqrt(marked_count * 1.0 / database_size));
    double b = (2.0 * iterations_num + 1.0) * a;
    double c = pow(sin(b), 2);
    return c;
}

int apply_Grover(int database_size, int* marked_indices, int marked_indices_num, int iterations_num) {
    // Check validity
    for (int i = 0; i < marked_indices_num; ++i) {
        assert(marked_indices[i] < database_size);
    }

    // Prepare the qubits to index all the elements in the database, and 1 extra for mark.
    int qubits_num = (int)ceil(log2(database_size)) + 1;
    Qureg* qureg = allocate_qureg(qubits_num);
    Qubit* marked_qubit = qureg->qubits[qubits_num - 1];
    Qubit** database_register = (Qubit**)pointer_memory_get(qubits_num - 1);
    for (int i = 0; i < qubits_num - 1; ++i) {
        database_register[i] = qureg->qubits[i];
    }
#if 0
    // State preparation (H) and Oracle
    apply_to_each(H.apply, database_register, qubits_num - 1);
    for (int i = 0; i < marked_indices_num; ++i) {
        CNOT.apply_on_int(qureg->qubits, qubits_num, marked_indices[i]);
    }

    for (int i = 0; i < iterations_num; ++i) {
        // Reflect marked state
        Z.apply(marked_qubit);

        // Adjoint state preparation (H) and Oracle
        for (int j = 0; j < marked_indices_num; ++j) {
            CNOT.apply_on_int_dagger(qureg->qubits, qubits_num, marked_indices[j]);
        }
        apply_to_each(H.apply_dagger, database_register, qubits_num - 1);
        printf("Adjoint state preparation (H) and Oracle\n");
        print_qureg(qureg);

        // Reflect zero state
        apply_to_each(X.apply, qureg->qubits, qubits_num);
        printf("Apply X to each qubits\n");
        print_qureg(qureg);
        CZ.apply(qureg->qubits, qubits_num);
        printf("Apply CZ to all qubits\n");
        print_qureg(qureg);
        apply_to_each(X.apply, qureg->qubits, qubits_num);
        printf("Reflect zero state\n");
        print_qureg(qureg);

        // State preparation (H) and Oracle
        apply_to_each(H.apply, database_register, qubits_num - 1);
        for (int j = 0; j < marked_indices_num; ++j) {
            CNOT.apply_on_int(qureg->qubits, qubits_num, marked_indices[j]);
        }
        printf("State preparation (H) and Oracle\n");
        print_qureg(qureg);
    }
#endif
#if 0
    X.apply(marked_qubit);
    apply_to_each(H.apply, qureg->qubits, qubits_num);
    for (int i = 0; i < iterations_num; ++i) {
        for (int i = 0; i < marked_indices_num; ++i) {
            CNOT.apply_on_int(qureg->qubits, qubits_num, marked_indices[i]);
        }
        apply_to_each(H.apply, database_register, qubits_num - 1);
        apply_to_each(X.apply, database_register, qubits_num - 1);
        CZ.apply(database_register, qubits_num - 1);
        apply_to_each_reverse(X.apply_dagger, database_register, qubits_num - 1);
        apply_to_each_reverse(H.apply_dagger, database_register, qubits_num - 1);
    }
#endif
    apply_X(marked_qubit);
    for (int i = 0; i < qubits_num; ++i) {
        apply_H(qureg->qubits[i]);
    }
    //print_qureg(qureg);
    for (int j = 0; j < iterations_num; ++j) {
        for (int i = 0; i < marked_indices_num; ++i) {
            apply_CNOT_on_int(qureg->qubits, qubits_num, marked_indices[i]);
        }
        //print_qureg(qureg);
        for (int i = 0; i < qubits_num - 1; ++i) {
            apply_H(database_register[i]);
        }
        //print_qureg(qureg);
        apply_P(database_register, qubits_num - 1);
        //print_qureg(qureg);
        for (int i = qubits_num - 2; i >= 0; --i) {
            apply_H(database_register[i]);
        }
        //print_qureg(qureg);
    }

    print_qureg(qureg);

    // Measurements
    for (int i = 0; i < qubits_num; ++i) {
        apply_PauliZ_M(qureg->qubits[i]);
    }
    int result = results_as_int(database_register, qubits_num - 1);

    pointer_memory_return(database_register, qubits_num - 1);
    free_qureg(qureg);

    return result;
}

int find_minimum(int* database, int database_size) {
    int min_index = rand() % database_size;
    int* marked_elements = int_memory_get(database_size - 1);

    while (1) {
        int marked_count = 0;
        for (int i = 0; i < database_size; ++i) {
            if (database[i] < database[min_index]) {
                marked_elements[marked_count] = i;
                ++marked_count;
            }
        }
        if (marked_count == 0) {
            return min_index;
        }

        int iterations_num = (int)ceil(PI / 4.0 * sqrt(database_size * 1.0 / marked_count));
        double successProbability = calculate_success_probability(marked_count, database_size, iterations_num);
        while(successProbability < 0.2) {
            iterations_num += 2;
            successProbability = calculate_success_probability(marked_count, database_size, iterations_num);
        }
        int trial_count = 0;
        
        while(1) {
            int index = apply_Grover(database_size, marked_elements, marked_count, iterations_num);
            if (database[index] < database[min_index]) {
                min_index = index;
                break;
            }
            assert(++trial_count < 1.0 / successProbability * 5);
        }
    }
    int_memory_return(marked_elements, database_size - 1);

    return min_index;
}
