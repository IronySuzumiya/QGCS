#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>

#include "Qubit.h"
#include "Gate.h"
#include "Util.h"
#include "Tensor.h"
#include "Algorithm.h"

typedef controlled_gate_apply keep_zero_phase_gate_apply;

typedef struct _keep_zero_phase_gate {
    keep_zero_phase_gate_apply apply;
} KeepZeroPhaseGate;

void apply_phase(Qubit** qubits, int qubits_num) {
    Qureg* qureg = qubits[0]->qureg;
    for (int i = 1; i < qubits_num; ++i) {
        assert(qubits[i]->qureg == qureg);
        assert(!qubits[i]->measured);
    }
    // As usual, check entanglement first
    bool all_product = true;
    for (int i = 0; i < qubits_num; ++i) {
        if (qubits[i]->entangled) {
            all_product = false;
            break;
        }
    }
    gsl_vector_complex* combined_inputs;
    int states_num;
    // Only used in the second situation
    int all_involved_qubits_num;
    int* order;
    int qupairs_num;
    int product_qubits_num;

    // If all qubits are not entangled
    if (all_product) {
        // Make up the combined vector and apply CNOT gate on it
        combined_inputs = Kronecker_product_vv(qubits[0]->state, qubits[1]->state);
        for (int i = 2; i < qubits_num; ++i) {
            gsl_vector_complex* temp = Kronecker_product_vv(combined_inputs, qubits[i]->state);
            gsl_vector_complex_free(combined_inputs);
            combined_inputs = temp;
        }
        states_num = combined_inputs->size;
        for (int i = 1; i < states_num; ++i) {
            gsl_complex temp = gsl_vector_complex_get(combined_inputs, i);
            gsl_vector_complex_set(combined_inputs, i, gsl_complex_negative(temp));
        }
    }
    else {
        //   Here maybe some shuffles applied to qubits, because some of them are in their own qupairs
        // in which there are some spaces between their ids or they have different order compared to
        // that in the inputs.
        // Work out how many independent qupairs there are, and store them
        Qupair** qupairs = (Qupair**)malloc(sizeof(Qupair*) * qubits_num);
        qupairs_num = 0;
        for (int i = 0; i < qubits_num; ++i) {
            if (qubits[i]->entangled) {
                bool exist = false;
                for (int j = 0; j < qupairs_num; ++j) {
                    if (qupairs[j] == qubits[i]->qupair) {
                        exist = true;
                        break;
                    }
                }
                if (!exist) {
                    qupairs[qupairs_num] = qubits[i]->qupair;
                    ++qupairs_num;
                }
            }
        }
        // Find out the order of qubits in the combined inputs
        all_involved_qubits_num = 0;
        product_qubits_num = 0;
        for (int i = 0; i < qupairs_num; ++i) {
            all_involved_qubits_num += qupairs[i]->qubits_num;
        }
        for (int i = 0; i < qubits_num; ++i) {
            if (!qubits[i]->entangled) {
                ++all_involved_qubits_num;
                ++product_qubits_num;
            }
        }
        order = (int*)malloc(sizeof(int) * all_involved_qubits_num);
        Qubit** product_qubits = (Qubit**)malloc(sizeof(Qubit*) * product_qubits_num);
        int index = 0;
        int product_qubits_index = 0;
        for (int i = 0; i < qupairs_num; ++i) {
            for (int j = 0; j < qupairs[i]->qubits_num; ++j) {
                order[index] = qupairs[i]->qubits_indices[j];
                ++index;
            }
        }
        for (int i = 0; i < qubits_num; ++i) {
            if (!qubits[i]->entangled) {
                order[index] = qubits[i]->index;
                ++index;
                product_qubits[product_qubits_index] = qubits[i];
                ++product_qubits_index;
            }
        }
        // Prepare the combined inputs
        if (qupairs_num > 1) {
            combined_inputs = Kronecker_product_vv(qupairs[0]->state, qupairs[1]->state);
            for (int i = 2; i < qupairs_num; ++i) {
                gsl_vector_complex* temp = Kronecker_product_vv(combined_inputs, qupairs[i]->state);
                gsl_vector_complex_free(combined_inputs);
                combined_inputs = temp;
            }
        }
        else {
            combined_inputs = gsl_vector_complex_calloc(qupairs[0]->states_num);
            gsl_vector_complex_memcpy(combined_inputs, qupairs[0]->state);
        }
        if (product_qubits_num > 0) {
            gsl_vector_complex* another_combined_inputs;
            if (product_qubits_num > 1) {
                another_combined_inputs = Kronecker_product_vv(product_qubits[0]->state, product_qubits[1]->state);
                for (int i = 2; i < qupairs_num; ++i) {
                    gsl_vector_complex* temp = Kronecker_product_vv(another_combined_inputs, product_qubits[i]->state);
                    gsl_vector_complex_free(another_combined_inputs);
                    another_combined_inputs = temp;
                }
            }
            else {
                another_combined_inputs = gsl_vector_complex_calloc(2);
                gsl_vector_complex_memcpy(another_combined_inputs, product_qubits[0]->state);
            }
            gsl_vector_complex* temp = Kronecker_product_vv(combined_inputs, another_combined_inputs);
            gsl_vector_complex_free(combined_inputs);
            gsl_vector_complex_free(another_combined_inputs);
            combined_inputs = temp;
        }
        // Shuffle the combined inputs to make the order right
        for (int i = 0; i < qubits_num; ++i) {
            bool success = false;
            for (int j = i; j < all_involved_qubits_num; ++j) {
                if (order[j] == qubits[i]->index) {
                    vector_complex_positions_swap(combined_inputs, all_involved_qubits_num, j, i);
                    int temp = order[j];
                    order[j] = order[i];
                    order[i] = temp;
                    success = true;
                }
            }
            assert(success);
        }
        states_num = combined_inputs->size;
        if (all_involved_qubits_num > qubits_num) {
            // There are outside qubits involved. Put them to the end of the order array and apply I on them.
            int ancillas_num = all_involved_qubits_num - qubits_num;
            int ancillas_states_num = (int)pow(2, ancillas_num);
            for (int i = ancillas_states_num; i < states_num; ++i) {
                gsl_complex temp = gsl_vector_complex_get(combined_inputs, i);
                gsl_vector_complex_set(combined_inputs, i, gsl_complex_negative(temp));
            }
        }
        else {
            for (int i = 1; i < states_num; ++i) {
                gsl_complex temp = gsl_vector_complex_get(combined_inputs, i);
                gsl_vector_complex_set(combined_inputs, i, gsl_complex_negative(temp));
            }
        }
        free(qupairs);
        free(product_qubits);
    }

    // After story ;-)

    int new_qupairs_num;
    // Prepare a new qupair list for the qureg
    if (all_product) {
        new_qupairs_num = qureg->qupairs_num - qubits_num + 1;
    }
    else {
        new_qupairs_num = qureg->qupairs_num - qupairs_num - product_qubits_num + 1;
    }
    Qupair** new_qupairs = (Qupair**)malloc(sizeof(Qupair*) * new_qupairs_num);
    // Mark those "clean" qupairs, which means they are not relevant to the inputs
    bool* dirty = (bool*)malloc(sizeof(bool) * qureg->qupairs_num);
    memset(dirty, 0, sizeof(bool) * qureg->qupairs_num);
    if (all_product) {
        for (int i = 0; i < qubits_num; ++i) {
            dirty[qubits[i]->qupair->index] = true;
        }
    }
    else {
        for (int i = 0; i < all_involved_qubits_num; ++i) {
            dirty[qureg->qubits[order[i]]->qupair->index] = true;
        }
    }
    // Delete "dirty qupairs and move "clean" qupairs to the new right place
    int new_qupair_index = 0;
    for (int i = 0; i < qureg->qupairs_num; ++i) {
        if (!dirty[i]) {
            new_qupairs[new_qupair_index] = qureg->qupairs[i];
            new_qupairs[new_qupair_index]->index = new_qupair_index;
            ++new_qupair_index;
        }
        else {
            free_qupair(qureg->qupairs[i]);
        }
    }
    free(dirty);
    // the new qupair which combines all the inputs
    Qupair* new_qupair = (Qupair*)malloc(sizeof(Qupair));
    new_qupairs[new_qupair_index] = new_qupair;
    new_qupair->index = new_qupair_index;
    new_qupair->qureg = qureg;
    if (all_product) {
        new_qupair->qubits_num = qubits_num;
        new_qupair->qubits_indices = (int*)malloc(sizeof(int) * qubits_num);
        for (int i = 0; i < qubits_num; ++i) {
            new_qupair->qubits_indices[i] = qubits[i]->index;
            qubits[i]->qupair = new_qupair;
            qubits[i]->entangled = true;
        }
    }
    else {
        new_qupair->qubits_num = all_involved_qubits_num;
        new_qupair->qubits_indices = (int*)malloc(sizeof(int) * all_involved_qubits_num);
        for (int i = 0; i < all_involved_qubits_num; ++i) {
            new_qupair->qubits_indices[i] = order[i];
            qureg->qubits[order[i]]->qupair = new_qupair;
            qureg->qubits[order[i]]->entangled = true;
        }
        free(order);
    }
    new_qupair->states_num = states_num;
    new_qupair->state = combined_inputs;
    // Only delete the qupair list and do not touch the content (Qupair*) in it
    free(qureg->qupairs);
    qureg->qupairs = new_qupairs;
    qureg->qupairs_num = new_qupairs_num;

    //   Consider the situation that no any input become entangled after this CNOT operation,
    // this design may look quite low-efficient, because I assume all of them become entangled
    // and delete their old qupairs and combine them directly. That's true, but actually this
    // is only the first version, or, say, a demo. So, for simplicity, the way I deal with this
    // problem is kind of acceptable.
    check_entanglement(new_qupair);
}

static KeepZeroPhaseGate phase = {
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
    apply_to_each(H.apply, qureg->qubits, qubits_num);
    for (int i = 0; i < iterations_num; ++i) {
        for (int i = 0; i < marked_indices_num; ++i) {
            CNOT.apply_on_int(qureg->qubits, qubits_num, marked_indices[i]);
        }
        //print_qureg(qureg);
        apply_to_each(H.apply, database_register, qubits_num - 1);
        //print_qureg(qureg);
        phase.apply(database_register, qubits_num - 1);
        //print_qureg(qureg);
        apply_to_each_reverse(H.apply_dagger, database_register, qubits_num - 1);
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
        // In a small-sized database, this parameter need to be tuned
        int iterations_num = (int)floor(PI / 4.0 * sqrt(database_size * 1.0 / marked_count));
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
