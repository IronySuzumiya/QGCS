#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <memory.h>

#include "Gate.h"
#include "Tensor.h"
#include "Util.h"
#include "Feature.h"
#include "Const.h"

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

static void apply_X_optimized(Qubit* qubit) {
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

static void apply_X(Qubit* qubit) {
#ifdef _QGCS_OPT
    apply_X_optimized(qubit);
#else
    apply_unitary_matrix(X.matrix, qubit);
#endif
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

void apply_controlled_gate_on_int(Qubit** qubits, int qubits_num, int controlled, do_with_controlled func, double param) {
    Qureg* qureg = qubits[0]->qureg;

    for (int i = 1; i < qubits_num; ++i) {
        assert(qubits[i]->qureg == qureg);
        assert(!qubits[i]->measured);
    }
    assert(controlled < pow(2, qubits_num - 1));

    bool all_separable = true;
    for (int i = 0; i < qubits_num; ++i) {
        if (qubits[i]->entangled) {
            all_separable = false;
            break;
        }
    }

    gsl_vector_complex* combined_inputs;
    int new_qupairs_num;
    int all_involved_qubits_num;
    // Only used in the second situation
    int* order;

    // If all qubits are not entangled
    if (all_separable) {
        // Make up the combined vector and apply CNOT gate on it
        combined_inputs = combine_qubits(qubits, qubits_num);

        all_involved_qubits_num = qubits_num;
        new_qupairs_num = qureg->qupairs_num - qubits_num + 1;

        func(combined_inputs, controlled, all_involved_qubits_num, qubits_num, param);
    }
    else {
        // Work out how many exclusive qupairs there are, and store them
        Qupair** qupairs = (Qupair**)malloc(sizeof(Qupair*) * qubits_num);
        int qupairs_num = get_exclusive_qupairs_from_qubits(qubits, qubits_num, qupairs);

        // Find out the order of qubits in the combined inputs
        int separable_qubits_num = get_separable_qubits_num(qubits, qubits_num);
        all_involved_qubits_num = separable_qubits_num + get_entangled_qubits_num(qupairs, qupairs_num);

        new_qupairs_num = qureg->qupairs_num - qupairs_num - separable_qubits_num + 1;

        order = (int*)malloc(sizeof(int) * all_involved_qubits_num);
        unfold_qupairs_and_qubits(order, qupairs, qupairs_num, qubits, qubits_num);

        Qubit** separable_qubits = (Qubit**)malloc(sizeof(Qubit*) * separable_qubits_num);
        get_separable_qubits(separable_qubits, qubits, qubits_num);

        // Prepare the combined inputs
        combined_inputs = combine_qupairs(qupairs, qupairs_num);
        if (separable_qubits_num > 0) {
            gsl_vector_complex* another_combined_inputs = combine_qubits(separable_qubits, separable_qubits_num);

            gsl_vector_complex* temp = Kronecker_product_vv(combined_inputs, another_combined_inputs);
            gsl_vector_complex_free(combined_inputs);
            gsl_vector_complex_free(another_combined_inputs);

            combined_inputs = temp;
        }

        free(qupairs);
        free(separable_qubits);

        // Shuffle the combined inputs to make the order right
        for (int i = 0; i < qubits_num; ++i) {
            for (int j = i; j < all_involved_qubits_num; ++j) {
                if (order[j] == qubits[i]->index) {
                    vector_complex_positions_swap(combined_inputs, all_involved_qubits_num, j, i);
                    int temp = order[j];
                    order[j] = order[i];
                    order[i] = temp;
                }
            }
        }

        func(combined_inputs, controlled, all_involved_qubits_num, qubits_num, param);
    }

    // After story ;-)

    // Mark those "clean" qupairs, which means they are not relevant to the inputs
    bool* dirty = (bool*)malloc(sizeof(bool) * qureg->qupairs_num);
    memset(dirty, 0, sizeof(bool) * qureg->qupairs_num);
    if (all_separable) {
        for (int i = 0; i < qubits_num; ++i) {
            dirty[qubits[i]->qupair->index] = true;
        }
    }
    else {
        for (int i = 0; i < all_involved_qubits_num; ++i) {
            dirty[qureg->qubits[order[i]]->qupair->index] = true;
        }
    }

    // Prepare a new qupair list for the qureg
    Qupair** new_qupairs = (Qupair**)malloc(sizeof(Qupair*) * new_qupairs_num);

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
    new_qupair->qubits_num = all_involved_qubits_num;
    new_qupair->qubits_indices = (int*)malloc(sizeof(int) * all_involved_qubits_num);
    new_qupair->states_num = combined_inputs->size;
    new_qupair->state = combined_inputs;

    if (all_separable) {
        for (int i = 0; i < qubits_num; ++i) {
            new_qupair->qubits_indices[i] = qubits[i]->index;
            qubits[i]->qupair = new_qupair;
            qubits[i]->entangled = true;
        }
    }
    else {
        for (int i = 0; i < all_involved_qubits_num; ++i) {
            new_qupair->qubits_indices[i] = order[i];
            qureg->qubits[order[i]]->qupair = new_qupair;
            qureg->qubits[order[i]]->entangled = true;
        }
        free(order);
    }

    // Only delete the qupair list and do not touch the content (Qupair*) in it
    free(qureg->qupairs);
    qureg->qupairs = new_qupairs;
    qureg->qupairs_num = new_qupairs_num;

    check_entanglement(new_qupair);
}

static void do_X_with_controlled(gsl_vector_complex* combined_inputs, int controlled, int all_involved_qubits_num, int qubits_num, double dummy) {
    assert(all_involved_qubits_num >= qubits_num);

    if (all_involved_qubits_num == qubits_num) {
        vector_complex_swap(combined_inputs, controlled * 2, controlled * 2 + 1);
    }
    else {
        int ancillas_num = all_involved_qubits_num - qubits_num;
        int ancillas_states_num = (int)pow(2, ancillas_num);
        for (int i = 0; i < ancillas_states_num; ++i) {
            vector_complex_swap(combined_inputs, controlled * 2 * ancillas_states_num + i, (controlled * 2 + 1) * ancillas_states_num + i);
        }
    }
}

static void apply_CNOT_on_int(Qubit** qubits, int qubits_num, int controlled) {
    apply_controlled_gate_on_int(qubits, qubits_num, controlled, do_X_with_controlled, 0.0);
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

static void do_R_with_controlled(gsl_vector_complex* combined_inputs, int controlled, int all_involved_qubits_num, int qubits_num, double phi) {
    assert(all_involved_qubits_num >= qubits_num);

    if (all_involved_qubits_num == qubits_num) {
        gsl_complex temp = gsl_complex_mul(gsl_vector_complex_get(combined_inputs, controlled * 2 + 1), gsl_complex_polar(1.0, phi));
        gsl_vector_complex_set(combined_inputs, controlled * 2 + 1, temp);
    }
    else {
        int ancillas_num = all_involved_qubits_num - qubits_num;
        int ancillas_states_num = (int)pow(2, ancillas_num);
        for (int i = 0; i < ancillas_states_num; ++i) {
            gsl_complex temp = gsl_vector_complex_get(combined_inputs, (controlled * 2 + 1) * ancillas_states_num + i);
            temp = gsl_complex_mul(temp, gsl_complex_polar(1.0, phi));
            gsl_vector_complex_set(combined_inputs, (controlled * 2 + 1) * ancillas_states_num + i, temp);
        }
    }
}

static void apply_CR_on_int(Qubit** qubits, int qubits_num, int controlled, double phi) {
    apply_controlled_gate_on_int(qubits, qubits_num, controlled, do_R_with_controlled, phi);
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

int apply_to_each(Gate gate, Qubit** qubits, int qubits_num) {
    for (int i = 0; i < qubits_num; ++i) {
        gate.apply(qubits[i]);
    }
    return 0;
}

int apply_to_each_reverse(Gate gate, Qubit** qubits, int qubits_num) {
    for (int i = qubits_num - 1; i >= 0; --i) {
        gate.apply(qubits[i]);
    }
    return 0;
}
