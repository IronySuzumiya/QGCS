#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>

#include "Algorithm.h"
#include "Qubit.h"
#include "Gate.h"
#include "Util.h"
#include "Tensor.h"
#include "Feature.h"
#include "Const.h"

typedef controlled_gate_apply multi_input_gate_apply;

typedef struct _multi_input_gate {
    multi_input_gate_apply apply;
} MultiInputGate;

static void do_phase_with_controlled(gsl_vector_complex* combined_inputs, int controlled, int all_involved_qubits_num, int qubits_num, double dummy) {
    assert(all_involved_qubits_num >= qubits_num);

    if (all_involved_qubits_num == qubits_num) {
        for (size_t i = 1; i < combined_inputs->size; ++i) {
            gsl_complex temp = gsl_vector_complex_get(combined_inputs, i);
            gsl_vector_complex_set(combined_inputs, i, gsl_complex_negative(temp));
        }
    }
    else {
        int ancillas_num = all_involved_qubits_num - qubits_num;
        int ancillas_states_num = (int)pow(2, ancillas_num);
        for (size_t i = ancillas_states_num; i < combined_inputs->size; ++i) {
            gsl_complex temp = gsl_vector_complex_get(combined_inputs, i);
            gsl_vector_complex_set(combined_inputs, i, gsl_complex_negative(temp));
        }
    }
}

void apply_phase(Qubit** qubits, int qubits_num) {
    apply_controlled_gate_on_int(qubits, qubits_num, 0, do_phase_with_controlled, 0.0);
}

static MultiInputGate phase = {
    apply_phase
};

int apply_Grover(int database_size, int* marked_indices, int marked_indices_num, int iterations_num) {
    // Check validity
    for (int i = 0; i < marked_indices_num; ++i) {
        assert(marked_indices[i] < database_size);
    }

    // Prepare the qubits to index all the elements in the database, and 1 extra for mark.
    int qubits_num = (int)ceil(log2(database_size)) + 1;
    Qureg* qureg = allocate_qureg(qubits_num);
    Qubit* marked_qubit = qureg->qubits[qubits_num - 1];
    Qubit** database_register = (Qubit**)malloc(sizeof(Qubit*) * (qubits_num - 1));
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
    X.apply(marked_qubit);
    apply_to_each(H, qureg->qubits, qubits_num);
    for (int i = 0; i < iterations_num; ++i) {
        for (int i = 0; i < marked_indices_num; ++i) {
            CNOT.apply_on_int(qureg->qubits, qubits_num, marked_indices[i]);
        }
        //print_qureg(qureg);
        apply_to_each(H, database_register, qubits_num - 1);
        //print_qureg(qureg);
        phase.apply(database_register, qubits_num - 1);
        //print_qureg(qureg);
        apply_to_each_reverse(H, database_register, qubits_num - 1);
        //print_qureg(qureg);
    }

    print_qureg(qureg);

    // Measurements
    for (int i = 0; i < qubits_num; ++i) {
        PauliZ_M.apply(qureg->qubits[i]);
    }
    int result = results_as_int(database_register, qubits_num - 1);

    free(database_register);

    return result;
}

int find_minimum(int* database, int database_size) {
    int min_index = rand() % database_size;
    int* marked_elements = (int*)malloc(sizeof(int) * (database_size - 1));

    while (true) {
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

        int iterations_num = (int)round(PI / 4.0 * sqrt(database_size * 1.0 / marked_count));
        double successProbability =
            pow(sin((2.0 * iterations_num + 1.0) * asin(sqrt(marked_count * 1.0 / database_size))), 2.0);
        int trial_count = 0;
        
        while(true) {
            int index = apply_Grover(database_size, marked_elements, marked_count, iterations_num);
            if (database[index] < database[min_index]) {
                min_index = index;
                break;
            }
            assert(++trial_count < 1.0 / successProbability * 6);
        }
    }
    free(marked_elements);

    return min_index;
}
