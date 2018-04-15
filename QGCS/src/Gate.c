#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <memory.h>

#include "Gate.h"
#include "Tensor.h"
#include "Util.h"

static int _check_entanglement(Qupair* qupair) {
    struct {
        enum { unassigned, normal, infinite } type;
        gsl_complex value;
    } rate = { 0 };
    gsl_complex temp_value = gsl_complex_rect(0, 0);
    int unentangled_index = -1;
    int states_num = qupair->states_num;

    // Qupair with single qubit means it's not entangled
    if (states_num == 2) {
        Qubit* left_qubit = qupair->qureg->qubits[qupair->qubits_indices[0]];
        left_qubit->entangled = false;
        gsl_vector_complex_memcpy(left_qubit->state, qupair->state);
        return 0;
    }

    for (int gap = states_num / 2, index = 0; gap >= 1; gap /= 2, ++index) {
        rate.type = unassigned;
        rate.value = gsl_complex_rect(0.0, 0.0);
        bool entangled = false;
        for (int times = 0; times < states_num / 2 / gap; times++) {
            for (int i = 0; i < gap; ++i) {
                int offset = times * 2 * gap;
                gsl_complex a = gsl_vector_complex_get(qupair->state, offset + i);
                gsl_complex b = gsl_vector_complex_get(qupair->state, offset + gap + i);
                bool a_zero = complex_is_zero(a);
                bool b_zero = complex_is_zero(b);
                if (rate.type == unassigned) {
                    if (a_zero && b_zero) {
                        //can be any
                    }
                    else if (b_zero) {
                        temp_value = a;
                        rate.type = infinite;
                    }
                    else {
                        if (a_zero) {
                            temp_value = b;
                        }
                        rate.type = normal;
                        rate.value = gsl_complex_div(a, b);
                    }
                }
                else {
                    if (a_zero && b_zero) {
                        //always equal
                    }
                    else if (b_zero) {
                        if (rate.type != infinite || !complex_equal(a, temp_value)) {
                            entangled = true;
                            break;
                        }
                    }
                    else {
                        if ((!complex_is_zero(temp_value) && a_zero && !complex_equal(b, temp_value))
                            || rate.type != normal || !complex_equal(rate.value, gsl_complex_div(a, b))) {
                            entangled = true;
                            break;
                        }
                    }
                }
            }
        }
        if (!entangled) {
            unentangled_index = index;
            goto unentangling;
        }
    }

unentangling:
    if (unentangled_index != -1) {
        assert(rate.type != unassigned);
        gsl_complex one = gsl_complex_rect(1.0, 0.0);
        gsl_complex zero = gsl_complex_rect(0.0, 0.0);
        gsl_complex alpha;
        gsl_complex beta;
        if (rate.type == normal) {
            gsl_complex r2_plus_1 = gsl_complex_add(gsl_complex_mul(rate.value, rate.value), one);
            gsl_complex sqrt_1_over_r2_plus_1 = gsl_complex_sqrt(gsl_complex_div(one, r2_plus_1));
            alpha = gsl_complex_mul(rate.value, sqrt_1_over_r2_plus_1);
            beta = sqrt_1_over_r2_plus_1;
        }
        else {
            alpha = one;
            beta = zero;
        }

        // Update the newly-unentangled qubit
        Qubit* unentangled_qubit = qupair->qureg->qubits[qupair->qubits_indices[unentangled_index]];
        unentangled_qubit->entangled = false;
        assert(!unentangled_qubit->measured);
        gsl_vector_complex_set(unentangled_qubit->state, 0, alpha);
        gsl_vector_complex_set(unentangled_qubit->state, 1, beta);

        // Update the qupair which this qubit belongs to
        Qupair* old_qupair = qupair;
        --old_qupair->qubits_num;
        old_qupair->states_num /= 2;
        // Update qubits indices
        // Avoid using realloc()
        int* new_qubits_indices = (int*)malloc(sizeof(int) * old_qupair->qubits_num);
        memcpy(new_qubits_indices, old_qupair->qubits_indices, sizeof(int) * unentangled_index);
        memcpy(&new_qubits_indices[unentangled_index],
            &old_qupair->qubits_indices[unentangled_index + 1], sizeof(int) * (old_qupair->qubits_num - unentangled_index));
        /*int new_qubits_index = 0;
        for (int i = 0; i < unentangled_index; ++i, ++new_qubits_index) {
        new_qubits_indices[new_qubits_index] = old_qupair->qubits_indices[i];
        }
        for (int i = unentangled_index + 1; i < old_qupair->qubits_num + 1; ++i, ++new_qubits_index) {
        new_qubits_indices[new_qubits_index] = old_qupair->qubits_indices[i];
        }*/
        free(old_qupair->qubits_indices);
        old_qupair->qubits_indices = new_qubits_indices;
        // Update qupair state
        int gap = states_num / (int)pow(2, unentangled_index + 1);
        int new_index = 0;
        gsl_vector_complex* new_state = gsl_vector_complex_calloc(states_num / 2);
        for (int times = 0; times < states_num / 2 / gap; times++) {
            for (int i = 0; i < gap; ++i) {
                int offset = times * 2 * gap;
                if (!complex_is_zero(alpha)) {
                    gsl_complex old_value = gsl_vector_complex_get(qupair->state, offset + i);
                    gsl_vector_complex_set(new_state, new_index, gsl_complex_div(old_value, alpha));
                }
                else {
                    gsl_complex old_value = gsl_vector_complex_get(qupair->state, offset + gap + i);
                    gsl_vector_complex_set(new_state, new_index, gsl_complex_div(old_value, beta));
                }
                ++new_index;
            }
        }
        gsl_vector_complex_free(old_qupair->state);
        old_qupair->state = new_state;

        // Create a new qupair with single qubit, and put it into the qureg
        Qureg* old_qureg = unentangled_qubit->qureg;
        ++old_qureg->qupairs_num;
        Qupair** new_qupairs = (Qupair**)malloc(sizeof(Qupair*) * old_qureg->qupairs_num);
        for (int i = 0; i < old_qureg->qupairs_num - 1; ++i) {
            new_qupairs[i] = old_qureg->qupairs[i];
        }
        new_qupairs[old_qureg->qupairs_num - 1] = (Qupair*)malloc(sizeof(Qupair));
        Qupair* new_qupair = new_qupairs[old_qureg->qupairs_num - 1];
        new_qupair->index = old_qureg->qupairs_num - 1;
        new_qupair->qureg = old_qureg;
        new_qupair->qubits_num = 1;
        new_qupair->qubits_indices = (int*)malloc(sizeof(int));
        new_qupair->qubits_indices[0] = unentangled_qubit->index;
        new_qupair->states_num = 2;
        new_qupair->state = gsl_vector_complex_calloc(2);
        gsl_vector_complex_set(new_qupair->state, 0, alpha);
        gsl_vector_complex_set(new_qupair->state, 1, beta);
        unentangled_qubit->qupair = new_qupair;
        // Only delete old qupair array, don't touch the elements within
        free(old_qureg->qupairs);
        old_qureg->qupairs = new_qupairs;

        // There maybe other unentangled qubits, so go on to check again
        return 1;
    }
    else {
        // All qubits in this qupair are entangled
        return 0;
    }
}

void check_entanglement(Qupair* qupair) {
    while (_check_entanglement(qupair));
}

static void apply_unitary_matrix(gsl_matrix_complex* matrix, Qubit* qubit) {
    assert(!qubit->measured);
    assert(matrix->size1 == matrix->size2 && matrix->size1 == 2);
    if (!qubit->entangled) {
        gsl_vector_complex* temp = matrix_vector_complex_mul(matrix, qubit->state);
        gsl_vector_complex_free(qubit->state);
        qubit->state = temp;
        gsl_vector_complex* temp2 = gsl_vector_complex_calloc(2);
        gsl_vector_complex_memcpy(temp2, temp);
        gsl_vector_complex_free(qubit->qupair->state);
        qubit->qupair->state = temp2;
    }
    else {
        int qubits_num = qubit->qupair->qubits_num;
        int states_num = qubit->qupair->states_num;
        int significance = qubit_index_in_qupair(qubit, qubit->qupair);
        int identity1_size = (int)pow(2, significance);
        int identity2_size = (int)pow(2, qubits_num - significance - 1);
        gsl_matrix_complex* final_matrix = gsl_matrix_complex_calloc(states_num, states_num);
        gsl_matrix_complex* identity1 = gsl_matrix_complex_calloc(identity1_size, identity1_size);
        gsl_matrix_complex* identity2 = gsl_matrix_complex_calloc(identity2_size, identity2_size);
        gsl_matrix_complex_set_identity(identity1);
        gsl_matrix_complex_set_identity(identity2);
        gsl_matrix_complex* middle_matrix = Kronecker_product_mm(identity1, matrix);
        final_matrix = Kronecker_product_mm(middle_matrix, identity2);
        gsl_vector_complex* result = matrix_vector_complex_mul(final_matrix, qubit->qupair->state);
        gsl_vector_complex_free(qubit->qupair->state);
        qubit->qupair->state = result;
        gsl_matrix_complex_free(final_matrix);
        gsl_matrix_complex_free(middle_matrix);
        gsl_matrix_complex_free(identity2);
        gsl_matrix_complex_free(identity1);

        check_entanglement(qubit->qupair);
    }
}

static void apply_H(Qubit* qubit) {
    apply_unitary_matrix(H.matrix, qubit);
}

static void apply_H_dagger(Qubit* qubit) {
    apply_H(qubit);
}

Gate H = {
    apply_H,
    apply_H_dagger,
    NULL,
    NULL
};

static void apply_X(Qubit* qubit) {
    assert(!qubit->measured);
    if (!qubit->entangled) {
        gsl_complex temp = gsl_vector_complex_get(qubit->state, 0);
        gsl_vector_complex_set(qubit->state, 0, gsl_vector_complex_get(qubit->state, 1));
        gsl_vector_complex_set(qubit->state, 1, temp);
    }
    else {
        int qubits_num = qubit->qupair->qubits_num;
        int states_num = qubit->qupair->states_num;
        int a = qubit_index_in_qupair(qubit, qubit->qupair);
        int b = qubits_num - 1 - a;
        //   0 I   0 0 ...
        //   I 0   0 0 ...
        //
        //   0 0   0 I ...
        //   0 0   I 0 ...
        //   ...   ... ...
        //
        // 2^(b+1) elements per cell
        // 2^a cells
        //
        // 2^b,   2^b+1,   2^b+2,   ... 2^(b+1) - 1, 0, 1, 2, ... 2^b-1,
        // 3*2^b, 3*2^b+1, 3*2^b+2, ...
        // 5*2^b, 5*2^b+1, ...
        int elements_per_cell = (int)pow(2, b + 1);
        int cells = (int)pow(2, a);
        int* new_order = (int*)malloc(sizeof(int) * states_num);
        int new_index = 0;
        for (int i = 0; i < cells; ++i) {
            for (int j = 0; j < elements_per_cell / 2; ++j) {
                new_order[new_index] = (1 + 2 * i) * elements_per_cell / 2 + j;
                ++new_index;
            }
            for (int j = 0; j < elements_per_cell / 2; ++j) {
                new_order[new_index] = (2 * i) * elements_per_cell / 2 + j;
                ++new_index;
            }
        }
        gsl_vector_complex* new_state = gsl_vector_complex_calloc(states_num);
        for (int i = 0; i < states_num; ++i) {
            gsl_vector_complex_set(new_state, i, gsl_vector_complex_get(qubit->qupair->state, new_order[i]));
        }
        free(new_order);
        gsl_vector_complex_free(qubit->qupair->state);
        qubit->qupair->state = new_state;

        check_entanglement(qubit->qupair);
    }
}

static void apply_X_dagger(Qubit* qubit) {
    apply_X(qubit);
}

Gate X = {
    apply_X,
    apply_X_dagger,
    NULL,
    NULL
};

static void apply_R(Qubit* qubit, double phi) {
    assert(!qubit->measured);
    if (!qubit->entangled) {
        gsl_complex temp = gsl_complex_mul(gsl_vector_complex_get(qubit->state, 1), gsl_complex_polar(1.0, phi));
        gsl_vector_complex_set(qubit->state, 1, temp);
    }
    else {
        int qubits_num = qubit->qupair->qubits_num;
        int states_num = qubit->qupair->states_num;
        int left = qubit_index_in_qupair(qubit, qubit->qupair);
        int left_id_size = (int)pow(2, left);
        int right = qubits_num - 1 - left;
        int right_id_size = (int)pow(2, right);
        //   I  0  0  0 ...
        //   0 -I  0  0 ...
        //
        //   0  0  I  0 ...
        //   0  0  0 -I ...
        //   ....  .... ...
        //
        //   2 * right_id_size elements per cell, left_id_size cells on aggregate
        int elements_per_cell = 2 * right_id_size;
        int cells = left_id_size;
        for (int i = 0; i < cells; ++i) {
            for (int j = 0; j < elements_per_cell / 2; ++j) {
                int index = i * elements_per_cell + elements_per_cell / 2 + j;
                gsl_complex temp = gsl_vector_complex_get(qubit->qupair->state, index);
                temp = gsl_complex_mul(temp, gsl_complex_polar(1.0, phi));
                gsl_vector_complex_set(qubit->qupair->state, index, temp);
            }
        }

        check_entanglement(qubit->qupair);
    }
}

static void apply_R_dagger(Qubit* qubit, double phi) {
    apply_R(qubit, -phi);
}

Gate1 R = {
    apply_R,
    apply_R_dagger
};

static void apply_Z(Qubit* qubit) {
    apply_R(qubit, PI);
}

static void apply_Z_dagger(Qubit* qubit) {
    apply_Z(qubit);
}

Gate Z = {
    apply_Z,
    apply_Z_dagger,
    NULL,
    NULL
};

static void apply_CNOT_on_int(Qubit** qubits, int qubits_num, int controlled) {
    Qureg* qureg = qubits[0]->qureg;
    for (int i = 1; i < qubits_num; ++i) {
        assert(qubits[i]->qureg == qureg);
        assert(!qubits[i]->measured);
    }
    assert(controlled < pow(2, qubits_num - 1));
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
        gsl_complex temp = gsl_vector_complex_get(combined_inputs, controlled * 2);
        gsl_vector_complex_set(combined_inputs, controlled * 2, gsl_vector_complex_get(combined_inputs, controlled * 2 + 1));
        gsl_vector_complex_set(combined_inputs, controlled * 2 + 1, temp);
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
            for (int i = 0; i < ancillas_states_num; ++i) {
                gsl_complex temp1 = gsl_vector_complex_get(combined_inputs, controlled * 2 * ancillas_states_num + i);
                gsl_complex temp2 = gsl_vector_complex_get(combined_inputs, (controlled * 2 + 1) * ancillas_states_num + i);
                gsl_vector_complex_set(combined_inputs, (controlled * 2 + 1) * ancillas_states_num + i, temp1);
                gsl_vector_complex_set(combined_inputs, controlled * 2 * ancillas_states_num + i, temp2);
            }
        }
        else {
            gsl_complex temp = gsl_vector_complex_get(combined_inputs, controlled * 2);
            gsl_vector_complex_set(combined_inputs, controlled * 2, gsl_vector_complex_get(combined_inputs, controlled * 2 + 1));
            gsl_vector_complex_set(combined_inputs, controlled * 2 + 1, temp);
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

static void apply_CNOT_on_int_dagger(Qubit** qubits, int qubits_num, int controlled) {
    apply_CNOT_on_int(qubits, qubits_num, controlled);
}

static void apply_CNOT(Qubit** qubits, int qubits_num) {
    apply_CNOT_on_int(qubits, qubits_num, (int)(pow(2, qubits_num - 1) - 1));
}

static void apply_CNOT_dagger(Qubit** qubits, int qubits_num) {
    apply_CNOT(qubits, qubits_num);
}

ControlledGate CNOT = {
    apply_CNOT,
    apply_CNOT_dagger,
    apply_CNOT_on_int,
    apply_CNOT_on_int_dagger
};

static void apply_CR_on_int(Qubit** qubits, int qubits_num, int controlled, double phi) {
    Qureg* qureg = qubits[0]->qureg;
    for (int i = 1; i < qubits_num; ++i) {
        assert(qubits[i]->qureg == qureg);
        assert(!qubits[i]->measured);
    }
    assert(controlled < pow(2, qubits_num - 1));
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
        gsl_complex temp = gsl_complex_mul(gsl_vector_complex_get(combined_inputs, controlled * 2 + 1), gsl_complex_polar(1.0, phi));
        gsl_vector_complex_set(combined_inputs, controlled * 2 + 1, temp);
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
            for (int i = 0; i < ancillas_states_num; ++i) {
                gsl_complex temp = gsl_vector_complex_get(combined_inputs, (controlled * 2 + 1) * ancillas_states_num + i);
                temp = gsl_complex_mul(temp, gsl_complex_polar(1.0, phi));
                gsl_vector_complex_set(combined_inputs, (controlled * 2 + 1) * ancillas_states_num + i, temp);
            }
        }
        else {
            gsl_complex temp = gsl_vector_complex_get(combined_inputs, controlled * 2 + 1);
            temp = gsl_complex_mul(temp, gsl_complex_polar(1.0, phi));
            gsl_vector_complex_set(combined_inputs, controlled * 2 + 1, temp);
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

static void apply_CR_on_int_dagger(Qubit** qubits, int qubits_num, int controlled, double phi) {
    apply_CR_on_int(qubits, qubits_num, controlled, -phi);
}

static void apply_CR(Qubit** qubits, int qubits_num, double phi) {
    apply_CR_on_int(qubits, qubits_num, (int)(pow(2, qubits_num - 1) - 1), phi);
}

static void apply_CR_dagger(Qubit** qubits, int qubits_num, double phi) {
    apply_CR(qubits, qubits_num, -phi);
}

ControlledGate1 CR = {
    apply_CR,
    apply_CR_dagger,
    apply_CR_on_int,
    apply_CR_on_int_dagger
};

static void apply_CZ_on_int(Qubit** qubits, int qubits_num, int controlled) {
    apply_CR_on_int(qubits, qubits_num, controlled, PI);
}

static void apply_CZ_on_int_dagger(Qubit** qubits, int qubits_num, int controlled) {
    apply_CZ_on_int(qubits, qubits_num, controlled);
}

static void apply_CZ(Qubit** qubits, int qubits_num) {
    apply_CR(qubits, qubits_num, PI);
}

static void apply_CZ_dagger(Qubit** qubits, int qubits_num) {
    apply_CZ(qubits, qubits_num);
}

ControlledGate CZ = {
    apply_CZ,
    apply_CZ_dagger,
    apply_CZ_on_int,
    apply_CZ_on_int_dagger
};

static void apply_PauliZ_M(Qubit* qubit) {
    if (qubit->measured) {
        // Repeated measurements do not change the result
        return;
    }
    if (!qubit->entangled) {
        double p_zero = gsl_complex_abs2(gsl_vector_complex_get(qubit->state, 0));
        double p_one = gsl_complex_abs2(gsl_vector_complex_get(qubit->state, 1));
        assert(double_equal(p_zero + p_one, 1.0));
        // For now
        qubit->value = double_is_zero(p_zero) ? One
            : double_is_zero(p_one) ? Zero
            : rand() / (RAND_MAX + 1.0) < p_zero ? Zero : One;
        qubit->measured = true;
    }
    else {
        int sig = qubit_index_in_qupair(qubit, qubit->qupair);
        assert(sig != -1);
        int* indices_where_qubit_is_zero = (int*)malloc(sizeof(int) * qubit->qupair->states_num / 2);
        int index_zero = 0;
        int* indices_where_qubit_is_one = (int*)malloc(sizeof(int) * qubit->qupair->states_num / 2);
        int index_one = 0;
        for (int i = 0; i < qubit->qupair->states_num; ++i) {
            if (i & (0x1 << (qubit->qupair->qubits_num - 1 - sig))) {
                indices_where_qubit_is_one[index_one] = i;
                ++index_one;
            }
            else {
                indices_where_qubit_is_zero[index_zero] = i;
                ++index_zero;
            }
        }
        assert(index_one == index_zero && index_one == qubit->qupair->states_num / 2);
        double p_zero = 0.0;
        double p_one = 0.0;
        for (int i = 0; i < qubit->qupair->states_num / 2; ++i) {
            p_zero += gsl_complex_abs2(gsl_vector_complex_get(qubit->qupair->state, indices_where_qubit_is_zero[i]));
            p_one += gsl_complex_abs2(gsl_vector_complex_get(qubit->qupair->state, indices_where_qubit_is_one[i]));
        }
        assert(double_equal(p_zero + p_one, 1.0));
        // For now
        qubit->value = double_is_zero(p_zero) ? One
            : double_is_zero(p_one) ? Zero
            : rand() / (RAND_MAX + 1.0) < p_zero ? Zero : One;
        qubit->measured = true;
        qubit->entangled = false;
        gsl_vector_complex* phi_prime = gsl_vector_complex_calloc(qubit->qupair->states_num / 2);
        if (qubit->value == Zero) {
            for (int i = 0; i < qubit->qupair->states_num / 2; ++i) {
                gsl_complex v = gsl_vector_complex_get(qubit->qupair->state, indices_where_qubit_is_zero[i]);
                v = gsl_complex_div(v, gsl_complex_rect(sqrt(p_zero), 0));
                gsl_vector_complex_set(phi_prime, i, v);
            }
        }
        else {
            for (int i = 0; i < qubit->qupair->states_num / 2; ++i) {
                gsl_complex v = gsl_vector_complex_get(qubit->qupair->state, indices_where_qubit_is_one[i]);
                v = gsl_complex_div(v, gsl_complex_rect(sqrt(p_one), 0));
                gsl_vector_complex_set(phi_prime, i, v);
            }
        }
        int* new_qubits_indices = (int*)malloc(sizeof(int) * (qubit->qupair->qubits_num - 1));
        for (int i = 0; i < sig; ++i) {
            new_qubits_indices[i] = qubit->qupair->qubits_indices[i];
        }
        for (int i = sig + 1; i < qubit->qupair->qubits_num; ++i) {
            new_qubits_indices[i - 1] = qubit->qupair->qubits_indices[i];
        }
        free(qubit->qupair->qubits_indices);
        qubit->qupair->qubits_indices = new_qubits_indices;
        --qubit->qupair->qubits_num;
        qubit->qupair->states_num /= 2;
        qubit->qupair->state = phi_prime;

        check_entanglement(qubit->qupair);
    }
}

Measurement PauliZ_M = {
    apply_PauliZ_M,
    NULL
};

int apply_to_each(gate_apply apply, Qubit** qubits, int qubits_num) {
    for (int i = 0; i < qubits_num; ++i) {
        apply(qubits[i]);
    }
    return 0;
}

int apply_to_each_reverse(gate_apply apply, Qubit** qubits, int qubits_num) {
    for (int i = qubits_num - 1; i >= 0; --i) {
        apply(qubits[i]);
    }
    return 0;
}

int gate_init(unsigned int seed) {
    H.matrix = gsl_matrix_complex_calloc(2, 2);
    gsl_matrix_complex_set_all(H.matrix, gsl_complex_rect(1.0 / sqrt(2), 0));
    gsl_matrix_complex_set(H.matrix, 1, 1, gsl_complex_rect(-1.0 / sqrt(2), 0));

    H.dagger_matrix = gsl_matrix_complex_calloc(2, 2);
    gsl_matrix_complex_memcpy(H.dagger_matrix, H.matrix);

    X.matrix = gsl_matrix_complex_calloc(2, 2);
    gsl_matrix_complex_set(X.matrix, 0, 1, gsl_complex_rect(1, 0));
    gsl_matrix_complex_set(X.matrix, 1, 0, gsl_complex_rect(1, 0));

    X.dagger_matrix = gsl_matrix_complex_calloc(2, 2);
    gsl_matrix_complex_memcpy(X.dagger_matrix, X.matrix);

    Z.matrix = gsl_matrix_complex_calloc(2, 2);
    gsl_matrix_complex_set(Z.matrix, 0, 0, gsl_complex_rect(1, 0));
    gsl_matrix_complex_set(Z.matrix, 1, 1, gsl_complex_rect(-1, 0));

    Z.dagger_matrix = gsl_matrix_complex_calloc(2, 2);
    gsl_matrix_complex_memcpy(Z.dagger_matrix, Z.matrix);

    PauliZ_M.operators = (gsl_matrix_complex**)malloc(sizeof(gsl_matrix_complex*) * 2);
    PauliZ_M.operators[0] = gsl_matrix_complex_calloc(2, 2);
    gsl_matrix_complex_set(PauliZ_M.operators[0], 0, 0, gsl_complex_rect(1.0, 0));
    PauliZ_M.operators[1] = gsl_matrix_complex_calloc(2, 2);
    gsl_matrix_complex_set(PauliZ_M.operators[1], 1, 1, gsl_complex_rect(1.0, 0));

    srand(seed);

    return 0;
}
