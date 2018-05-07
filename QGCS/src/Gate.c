#include <math.h>
#include <assert.h>

// This should be deleted in the future
#include <stdlib.h>

#include "Gate.h"
#include "Memory.h"
#include "Tensor.h"
#include "Util.h"
#include "Feature.h"
#include "Const.h"

static void apply_unitary_matrix(Matrix matrix, Qubit* qubit) {
    assert(!qubit->measured);
    assert(matrix.size1 == matrix.size2 && matrix.size1 == 2);
    if (!qubit->entangled) {
        Ket temp = matrix_mul_ket(matrix, qubit->state);
        ket_free(qubit->state);
        qubit->state = temp;
        Ket temp2 = ket_from_ket(temp);
        ket_free(qubit->qupair->state);
        qubit->qupair->state = temp2;
    }
    else {
        int qubits_num = qubit->qupair->qubits_num;
        int states_num = qubit->qupair->states_num;
        int significance = qubit_index_in_qupair(qubit, qubit->qupair);
        int identity1_size = (int)pow(2, significance);
        int identity2_size = (int)pow(2, qubits_num - significance - 1);
        Matrix final_matrix = matrix_calloc(states_num, states_num);
        Matrix identity1 = matrix_calloc(identity1_size, identity1_size);
        Matrix identity2 = matrix_calloc(identity2_size, identity2_size);
        matrix_set_identity(identity1);
        matrix_set_identity(identity2);
        Matrix middle_matrix = matrix_Kronecker_product(identity1, matrix);
        final_matrix = matrix_Kronecker_product(middle_matrix, identity2);
        Ket result = matrix_mul_ket(final_matrix, qubit->qupair->state);
        ket_free(qubit->qupair->state);
        qubit->qupair->state = result;
        matrix_free(final_matrix);
        matrix_free(middle_matrix);
        matrix_free(identity2);
        matrix_free(identity1);

        check_entanglement(qubit->qupair);
    }
}

void apply_H(Qubit* qubit) {
    apply_unitary_matrix(H, qubit);
}

void apply_H_dagger(Qubit* qubit) {
    apply_H(qubit);
}

static Complex H_block[4] = { { 0.70710678, 0 }, { 0.70710678, 0 }, { 0.70710678, 0 }, { -0.70710678, 0 } };

Matrix H = {
    H_block, 2, 2
};

void apply_X(Qubit* qubit) {
    assert(!qubit->measured);
    if (!qubit->entangled) {
        Complex temp = ket_get(qubit->state, 0);
        ket_set(qubit->state, 0, ket_get(qubit->state, 1));
        ket_set(qubit->state, 1, temp);
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
        int* new_order = int_memory_get(states_num);
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
        Ket new_state = ket_calloc(states_num);
        for (int i = 0; i < states_num; ++i) {
            ket_set(new_state, i, ket_get(qubit->qupair->state, new_order[i]));
        }
        int_memory_return(new_order, states_num);
        ket_free(qubit->qupair->state);
        qubit->qupair->state = new_state;

        check_entanglement(qubit->qupair);
    }
}

void apply_X_dagger(Qubit* qubit) {
    apply_X(qubit);
}

void apply_R(Qubit* qubit, double phi) {
    assert(!qubit->measured);
    if (!qubit->entangled) {
        Complex temp = complex_mul(ket_get(qubit->state, 1), complex_polar(1.0f, phi));
        ket_set(qubit->state, 1, temp);
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
                Complex temp = ket_get(qubit->qupair->state, index);
                temp = complex_mul(temp, complex_polar(1.0f, phi));
                ket_set(qubit->qupair->state, index, temp);
            }
        }

        check_entanglement(qubit->qupair);
    }
}

void apply_R_dagger(Qubit* qubit, double phi) {
    apply_R(qubit, -phi);
}

void apply_Z(Qubit* qubit) {
    apply_R(qubit, PI);
}

void apply_Z_dagger(Qubit* qubit) {
    apply_Z(qubit);
}

static int check_all_separable(Qubit** qubits, int qubits_num) {
    Qureg* qureg = qubits[0]->qureg;
    for (int i = 1; i < qubits_num; ++i) {
        assert(qubits[i]->qureg == qureg);
        assert(!qubits[i]->measured);
    }
    int all_separable = 1;
    for (int i = 0; i < qubits_num; ++i) {
        if (qubits[i]->entangled) {
            all_separable = 0;
            break;
        }
    }
    return all_separable;
}

static int* get_dirty_list_of_separable_qubits(Qubit** qubits, int qubits_num) {
    Qureg* qureg = qubits[0]->qureg;
    int* dirty = int_memory_get(qureg->qupairs_num);
    for (int i = 0; i < qureg->qupairs_num; ++i) {
        dirty[i] = 0;
    }
    for (int i = 0; i < qubits_num; ++i) {
        dirty[qubits[i]->qupair->index] = 1;
    }
    return dirty;
}

static int* get_dirty_list_of_mixed_qubits(Qubit** qubits, int* order, int all_involved_qubits_num) {
    Qureg* qureg = qubits[0]->qureg;
    int* dirty = int_memory_get(qureg->qupairs_num);
    for (int i = 0; i < qureg->qupairs_num; ++i) {
        dirty[i] = 0;
    }
    for (int i = 0; i < all_involved_qubits_num; ++i) {
        dirty[qureg->qubits[order[i]]->qupair->index] = 1;
    }
    return dirty;
}

static void prepare_new_qupair_list(Qureg* qureg, int* dirty, int new_qupairs_num) {
    Qupair** new_qupairs = (Qupair**)pointer_memory_get(new_qupairs_num);
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
    int_memory_return(dirty, qureg->qupairs_num);
    pointer_memory_return(qureg->qupairs, qureg->qupairs_num);
    qureg->qupairs = new_qupairs;
    qureg->qupairs_num = new_qupairs_num;
}

static void set_new_qupair_profile(Qupair* new_qupair, Qureg* qureg, int qubits_num, Ket combined_inputs) {
    qureg->qupairs[qureg->qupairs_num - 1] = new_qupair;
    new_qupair->index = qureg->qupairs_num - 1;
    new_qupair->qureg = qureg;
    new_qupair->qubits_num = qubits_num;
    new_qupair->qubits_indices = int_memory_get(qubits_num);
    new_qupair->states_num = combined_inputs.size;
    new_qupair->state = combined_inputs;
}

static void set_new_qupair_qubits_indices_of_separable_qubits(Qupair* new_qupair, Qubit** qubits, int qubits_num) {
    for (int i = 0; i < qubits_num; ++i) {
        new_qupair->qubits_indices[i] = qubits[i]->index;
        qubits[i]->qupair = new_qupair;
        qubits[i]->entangled = 1;
    }
}

static void set_new_qupair_qubits_indices_of_mixed_qubits(Qupair* new_qupair, Qureg* qureg, int* order, int all_involved_qubits_num) {
    for (int i = 0; i < all_involved_qubits_num; ++i) {
        new_qupair->qubits_indices[i] = order[i];
        qureg->qubits[order[i]]->qupair = new_qupair;
        qureg->qubits[order[i]]->entangled = 1;
    }
    int_memory_return(order, all_involved_qubits_num);
}

static void set_new_qupair_of_separable_qubits(Qubit** qubits, int qubits_num, Ket combined_inputs) {
    Qupair* new_qupair = (Qupair*)qupair_memory_get(1);
    set_new_qupair_profile(new_qupair, qubits[0]->qureg, qubits_num, combined_inputs);
    set_new_qupair_qubits_indices_of_separable_qubits(new_qupair, qubits, qubits_num);
}

static void set_new_qupair_of_mixed_qubits(Qubit** qubits, int* order, int all_involved_qubits_num, Ket combined_inputs) {
    Qupair* new_qupair = (Qupair*)qupair_memory_get(1);
    set_new_qupair_profile(new_qupair, qubits[0]->qureg, all_involved_qubits_num, combined_inputs);
    set_new_qupair_qubits_indices_of_mixed_qubits(new_qupair, qubits[0]->qureg, order, all_involved_qubits_num);
}

static void set_new_qupair_list_of_separable_qubits(Qubit** qubits, int qubits_num, int new_qupairs_num, Ket combined_inputs) {
    Qureg* qureg = qubits[0]->qureg;
    int* dirty = get_dirty_list_of_separable_qubits(qubits, qubits_num);
    prepare_new_qupair_list(qureg, dirty, new_qupairs_num);
    set_new_qupair_of_separable_qubits(qubits, qubits_num, combined_inputs);
}

static void set_new_qupair_list_of_mixed_qubits(Qubit** qubits, int* order, int all_involved_qubits_num, int new_qupairs_num, Ket combined_inputs) {
    Qureg* qureg = qubits[0]->qureg;
    int* dirty = get_dirty_list_of_mixed_qubits(qubits, order, all_involved_qubits_num);
    prepare_new_qupair_list(qureg, dirty, new_qupairs_num);
    set_new_qupair_of_mixed_qubits(qubits, order, all_involved_qubits_num, combined_inputs);
}

static void update_qureg_after_applying_gate_on_separable_qubits(Qubit** qubits, int qubits_num, Ket combined_inputs) {
    Qureg* qureg = qubits[0]->qureg;
    set_new_qupair_list_of_separable_qubits(qubits, qubits_num, qureg->qupairs_num - qubits_num + 1, combined_inputs);
    check_entanglement(qureg->qupairs[qureg->qupairs_num - 1]);
}

static int* shuffle_qubit_list(int all_involved_qubits_num, Qupair** qupairs, int qupairs_num, Qubit** qubits, int qubits_num, Ket combined_inputs) {
    int *order = int_memory_get(all_involved_qubits_num);
    unfold_qupairs_and_qubits(order, qupairs, qupairs_num, qubits, qubits_num);
    for (int i = 0; i < qubits_num; ++i) {
        for (int j = i; j < all_involved_qubits_num; ++j) {
            if (order[j] == qubits[i]->index) {
                ket_positions_swap(combined_inputs, all_involved_qubits_num, j, i);
                int temp = order[j];
                order[j] = order[i];
                order[i] = temp;
            }
        }
    }
    return order;
}

static Ket combine_mixed_qubits(Qubit** qubits, int qubits_num, int** order, int* all_involved_qubits_num, int* new_qupairs_num) {
    Qupair** qupairs = (Qupair**)pointer_memory_get(qubits_num);
    int qupairs_num = get_exclusive_qupairs_from_qubits(qubits, qubits_num, qupairs);
    int separable_qubits_num = get_separable_qubits_num(qubits, qubits_num);
    *new_qupairs_num = qubits[0]->qureg->qupairs_num - qupairs_num - separable_qubits_num + 1;
    *all_involved_qubits_num = separable_qubits_num + get_entangled_qubits_num(qupairs, qupairs_num);
    Ket combined_inputs = combine_qupairs(qupairs, qupairs_num);

    if (separable_qubits_num > 0) {
        Qubit** separable_qubits = (Qubit**)pointer_memory_get(separable_qubits_num);
        get_separable_qubits(separable_qubits, qubits, qubits_num);
        Ket another_combined_inputs = combine_qubits(separable_qubits, separable_qubits_num);
        Ket temp = ket_Kronecker_product(combined_inputs, another_combined_inputs);
        ket_free(combined_inputs);
        ket_free(another_combined_inputs);
        combined_inputs = temp;
        pointer_memory_return(separable_qubits, separable_qubits_num);
    }

    *order = shuffle_qubit_list(*all_involved_qubits_num, qupairs, qupairs_num, qubits, qubits_num, combined_inputs);

    pointer_memory_return(qupairs, qubits_num);

    return combined_inputs;
}

static void update_qureg_after_applying_gate_on_mixed_qubits(Qubit** qubits, int* order, int all_involved_qubits_num, int new_qupairs_num, Ket combined_inputs) {
    Qureg* qureg = qubits[0]->qureg;
    set_new_qupair_list_of_mixed_qubits(qubits, order, all_involved_qubits_num, new_qupairs_num, combined_inputs);
    check_entanglement(qureg->qupairs[qureg->qupairs_num - 1]);
}

static void do_X(Ket combined_inputs, int controlled, int all_involved_qubits_num, int qubits_num) {
    assert(all_involved_qubits_num >= qubits_num);

    if (all_involved_qubits_num == qubits_num) {
        ket_swap(combined_inputs, controlled * 2, controlled * 2 + 1);
    }
    else {
        int ancillas_num = all_involved_qubits_num - qubits_num;
        int ancillas_states_num = (int)pow(2, ancillas_num);
        for (int i = 0; i < ancillas_states_num; ++i) {
            ket_swap(combined_inputs, controlled * 2 * ancillas_states_num + i, (controlled * 2 + 1) * ancillas_states_num + i);
        }
    }
}

void apply_CNOT_on_int(Qubit** qubits, int qubits_num, int controlled) {
    assert(controlled < pow(2, qubits_num - 1));
    Ket combined_inputs;
    int all_involved_qubits_num;
    int new_qupairs_num;
    int* order = NULL;

    if (check_all_separable(qubits, qubits_num)) {
        combined_inputs = combine_qubits(qubits, qubits_num);
        do_X(combined_inputs, controlled, qubits_num, qubits_num);
        update_qureg_after_applying_gate_on_separable_qubits(qubits, qubits_num, combined_inputs);
    }
    else {
        combined_inputs = combine_mixed_qubits(qubits, qubits_num, &order, &all_involved_qubits_num, &new_qupairs_num);
        do_X(combined_inputs, controlled, all_involved_qubits_num, qubits_num);
        update_qureg_after_applying_gate_on_mixed_qubits(qubits, order, all_involved_qubits_num, new_qupairs_num, combined_inputs);
    }
}

void apply_CNOT_dagger_on_int(Qubit** qubits, int qubits_num, int controlled) {
    apply_CNOT_on_int(qubits, qubits_num, controlled);
}

void apply_CNOT(Qubit** qubits, int qubits_num) {
    apply_CNOT_on_int(qubits, qubits_num, (int)(pow(2, qubits_num - 1) - 1));
}

void apply_CNOT_dagger(Qubit** qubits, int qubits_num) {
    apply_CNOT(qubits, qubits_num);
}

static void do_R(Ket combined_inputs, int controlled, int all_involved_qubits_num, int qubits_num, double phi) {
    assert(all_involved_qubits_num >= qubits_num);

    if (all_involved_qubits_num == qubits_num) {
        Complex temp = complex_mul(ket_get(combined_inputs, controlled * 2 + 1), complex_polar(1.0f, phi));
        ket_set(combined_inputs, controlled * 2 + 1, temp);
    }
    else {
        int ancillas_num = all_involved_qubits_num - qubits_num;
        int ancillas_states_num = (int)pow(2, ancillas_num);
        for (int i = 0; i < ancillas_states_num; ++i) {
            Complex temp = ket_get(combined_inputs, (controlled * 2 + 1) * ancillas_states_num + i);
            temp = complex_mul(temp, complex_polar(1.0f, phi));
            ket_set(combined_inputs, (controlled * 2 + 1) * ancillas_states_num + i, temp);
        }
    }
}

void apply_CR_on_int(Qubit** qubits, int qubits_num, int controlled, double phi) {
    assert(controlled < pow(2, qubits_num - 1));
    Ket combined_inputs;
    int all_involved_qubits_num;
    int new_qupairs_num;
    int* order = NULL;

    if (check_all_separable(qubits, qubits_num)) {
        combined_inputs = combine_qubits(qubits, qubits_num);
        do_R(combined_inputs, controlled, qubits_num, qubits_num, phi);
        update_qureg_after_applying_gate_on_separable_qubits(qubits, qubits_num, combined_inputs);
    }
    else {
        combined_inputs = combine_mixed_qubits(qubits, qubits_num, &order, &all_involved_qubits_num, &new_qupairs_num);
        do_R(combined_inputs, controlled, all_involved_qubits_num, qubits_num, phi);
        update_qureg_after_applying_gate_on_mixed_qubits(qubits, order, all_involved_qubits_num, new_qupairs_num, combined_inputs);
    }
}

void apply_CR_dagger_on_int(Qubit** qubits, int qubits_num, int controlled, double phi) {
    apply_CR_on_int(qubits, qubits_num, controlled, -phi);
}

void apply_CR(Qubit** qubits, int qubits_num, double phi) {
    apply_CR_on_int(qubits, qubits_num, (int)(pow(2, qubits_num - 1) - 1), phi);
}

void apply_CR_dagger(Qubit** qubits, int qubits_num, double phi) {
    apply_CR(qubits, qubits_num, -phi);
}

static void do_Z(Ket combined_inputs, int controlled, int all_involved_qubits_num, int qubits_num) {
    assert(all_involved_qubits_num >= qubits_num);

    if (all_involved_qubits_num == qubits_num) {
        Complex temp = complex_neg(ket_get(combined_inputs, controlled * 2 + 1));
        ket_set(combined_inputs, controlled * 2 + 1, temp);
    }
    else {
        int ancillas_num = all_involved_qubits_num - qubits_num;
        int ancillas_states_num = (int)pow(2, ancillas_num);
        for (int i = 0; i < ancillas_states_num; ++i) {
            Complex temp = complex_neg(ket_get(combined_inputs, (controlled * 2 + 1) * ancillas_states_num + i));
            ket_set(combined_inputs, (controlled * 2 + 1) * ancillas_states_num + i, temp);
        }
    }
}

void apply_CZ_on_int(Qubit** qubits, int qubits_num, int controlled) {
    assert(controlled < pow(2, qubits_num - 1));
    Ket combined_inputs;
    int all_involved_qubits_num;
    int new_qupairs_num;
    int* order = NULL;

    if (check_all_separable(qubits, qubits_num)) {
        combined_inputs = combine_qubits(qubits, qubits_num);
        do_Z(combined_inputs, controlled, qubits_num, qubits_num);
        update_qureg_after_applying_gate_on_separable_qubits(qubits, qubits_num, combined_inputs);
    }
    else {
        combined_inputs = combine_mixed_qubits(qubits, qubits_num, &order, &all_involved_qubits_num, &new_qupairs_num);
        do_Z(combined_inputs, controlled, all_involved_qubits_num, qubits_num);
        update_qureg_after_applying_gate_on_mixed_qubits(qubits, order, all_involved_qubits_num, new_qupairs_num, combined_inputs);
    }
}

void apply_CZ_dagger_on_int(Qubit** qubits, int qubits_num, int controlled) {
    apply_CZ_on_int(qubits, qubits_num, controlled);
}

void apply_CZ(Qubit** qubits, int qubits_num) {
    apply_CZ_on_int(qubits, qubits_num, (int)(pow(2, qubits_num - 1) - 1));
}

void apply_CZ_dagger(Qubit** qubits, int qubits_num) {
    apply_CZ(qubits, qubits_num);
}

// Technically, this ought to be implemented on PC terminal
void apply_PauliZ_M(Qubit* qubit) {
    if (qubit->measured) {
        // Repeated measurements do not change the result
        return;
    }
    if (!qubit->entangled) {
        double p_zero = complex_norm(ket_get(qubit->state, 0));
        double p_one = complex_norm(ket_get(qubit->state, 1));
        assert(double_equal(p_zero + p_one, 1.0));
        
        qubit->value = double_is_zero(p_zero) ? One
            : double_is_zero(p_one) ? Zero
            : rand() / (RAND_MAX + 1.0) < p_zero ? Zero : One;
        qubit->measured = 1;
    }
    else {
        int sig = qubit_index_in_qupair(qubit, qubit->qupair);
        assert(sig != -1);
        int* indices_where_qubit_is_zero = int_memory_get(qubit->qupair->states_num / 2);
        int index_zero = 0;
        int* indices_where_qubit_is_one = int_memory_get(qubit->qupair->states_num / 2);
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
            p_zero += complex_norm(ket_get(qubit->qupair->state, indices_where_qubit_is_zero[i]));
            p_one += complex_norm(ket_get(qubit->qupair->state, indices_where_qubit_is_one[i]));
        }
        assert(double_equal(p_zero + p_one, 1.0));
        
        qubit->value = double_is_zero(p_zero) ? One
            : double_is_zero(p_one) ? Zero
            : rand() / (RAND_MAX + 1.0) < p_zero ? Zero : One;
        qubit->measured = 1;
        qubit->entangled = 0;

        Ket phi_prime = ket_calloc(qubit->qupair->states_num / 2);
        if (qubit->value == Zero) {
            for (int i = 0; i < qubit->qupair->states_num / 2; ++i) {
                Complex v = ket_get(qubit->qupair->state, indices_where_qubit_is_zero[i]);
                v = complex_div(v, complex_rect(sqrt(p_zero), 0));
                ket_set(phi_prime, i, v);
            }
        }
        else {
            for (int i = 0; i < qubit->qupair->states_num / 2; ++i) {
                Complex v = ket_get(qubit->qupair->state, indices_where_qubit_is_one[i]);
                v = complex_div(v, complex_rect(sqrt(p_one), 0));
                ket_set(phi_prime, i, v);
            }
        }
        int_memory_return(indices_where_qubit_is_zero, qubit->qupair->states_num / 2);
        int_memory_return(indices_where_qubit_is_one, qubit->qupair->states_num / 2);

        int* new_qubits_indices = int_memory_get(qubit->qupair->qubits_num - 1);
        for (int i = 0; i < sig; ++i) {
            new_qubits_indices[i] = qubit->qupair->qubits_indices[i];
        }
        for (int i = sig + 1; i < qubit->qupair->qubits_num; ++i) {
            new_qubits_indices[i - 1] = qubit->qupair->qubits_indices[i];
        }
        int_memory_return(qubit->qupair->qubits_indices, qubit->qupair->qubits_num);
        qubit->qupair->qubits_indices = new_qubits_indices;

        --qubit->qupair->qubits_num;
        qubit->qupair->states_num /= 2;
        qubit->qupair->state = phi_prime;

        check_entanglement(qubit->qupair);
    }
}

static void do_P(Ket combined_inputs, int all_involved_qubits_num, int qubits_num) {
    assert(all_involved_qubits_num >= qubits_num);

    if (all_involved_qubits_num == qubits_num) {
        for (int i = 1; i < combined_inputs.size; ++i) {
            Complex temp = ket_get(combined_inputs, i);
            ket_set(combined_inputs, i, complex_neg(temp));
        }
    }
    else {
        int ancillas_num = all_involved_qubits_num - qubits_num;
        int ancillas_states_num = (int)pow(2, ancillas_num);
        for (int i = ancillas_states_num; i < combined_inputs.size; ++i) {
            Complex temp = ket_get(combined_inputs, i);
            ket_set(combined_inputs, i, complex_neg(temp));
        }
    }
}

void apply_P(Qubit** qubits, int qubits_num) {
    Ket combined_inputs;
    int all_involved_qubits_num;
    int new_qupairs_num;
    int* order = NULL;

    if (check_all_separable(qubits, qubits_num)) {
        combined_inputs = combine_qubits(qubits, qubits_num);
        do_P(combined_inputs, qubits_num, qubits_num);
        update_qureg_after_applying_gate_on_separable_qubits(qubits, qubits_num, combined_inputs);
    }
    else {
        combined_inputs = combine_mixed_qubits(qubits, qubits_num, &order, &all_involved_qubits_num, &new_qupairs_num);
        do_P(combined_inputs, all_involved_qubits_num, qubits_num);
        update_qureg_after_applying_gate_on_mixed_qubits(qubits, order, all_involved_qubits_num, new_qupairs_num, combined_inputs);
    }
}

void apply_P_dagger(Qubit** qubits, int qubits_num) {
    apply_P(qubits, qubits_num);
}
