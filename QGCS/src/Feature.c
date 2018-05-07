#include <assert.h>
#include <math.h>

#include "Feature.h"
#include "Util.h"
#include "Tensor.h"
#include "Memory.h"

static int _check_entanglement(Qupair* qupair, Complex* ratio, int* is_normal) {
    enum { unassigned, normal, infinite } type;
    Complex temp_value = complex_rect(0, 0);
    int separable_index = -1;
    int states_num = qupair->states_num;

    for (int gap = states_num / 2, index = 0; gap >= 1; gap /= 2, ++index) {
        type = unassigned;
        *ratio = complex_rect(0.0, 0.0);
        int entangled = 0;
        for (int times = 0; times < states_num / 2 / gap; times++) {
            for (int i = 0; i < gap; ++i) {
                int offset = times * 2 * gap;
                Complex a = ket_get(qupair->state, offset + i);
                Complex b = ket_get(qupair->state, offset + gap + i);
                int a_zero = complex_is_zero(a);
                int b_zero = complex_is_zero(b);
                if (type == unassigned) {
                    if (a_zero && b_zero) {
                        //can be any
                    }
                    else if (b_zero) {
                        temp_value = a;
                        type = infinite;
                    }
                    else {
                        if (a_zero) {
                            temp_value = b;
                        }
                        type = normal;
                        *ratio = complex_div(a, b);
                    }
                }
                else {
                    if (a_zero && b_zero) {
                        //always equal
                    }
                    else if (b_zero) {
                        if (type != infinite || !complex_equal(a, temp_value)) {
                            entangled = 1;
                            break;
                        }
                    }
                    else {
                        if ((!complex_is_zero(temp_value) && a_zero && !complex_equal(b, temp_value))
                            || type != normal || !complex_equal(*ratio, complex_div(a, b))) {
                            entangled = 1;
                            break;
                        }
                    }
                }
            }
        }
        if (!entangled) {
            separable_index = index;
            assert(type != unassigned);
            *is_normal = type == normal;
            break;
        }
    }

    return separable_index;
}

static int calculate_probamp(Complex ratio, int is_normal, Complex* alpha, Complex* beta) {
    Complex one = complex_rect(1.0, 0.0);
    Complex zero = complex_rect(0.0, 0.0);
    if (is_normal) {
        Complex r2_plus_1 = complex_add(complex_mul(ratio, ratio), one);
        Complex sqrt_1_over_r2_plus_1 = complex_sqrt(complex_div(one, r2_plus_1));
        *alpha = complex_mul(ratio, sqrt_1_over_r2_plus_1);
        *beta = sqrt_1_over_r2_plus_1;
    }
    else {
        *alpha = one;
        *beta = zero;
    }

    return 0;
}

static int detach_qubit_from_qupair(Qupair* qupair, int index, Complex alpha, Complex beta) {
    int* new_qubits_indices = int_memory_get(qupair->qubits_num - 1);
    for (int i = 0; i < index; ++i) {
        new_qubits_indices[i] = qupair->qubits_indices[i];
    }
    for (int i = index; i < qupair->qubits_num - 1; ++i) {
        new_qubits_indices[i] = qupair->qubits_indices[i + 1];
    }
    int_memory_return(qupair->qubits_indices, qupair->qubits_num);
    --qupair->qubits_num;
    qupair->states_num /= 2;
    qupair->qubits_indices = new_qubits_indices;

    int gap = qupair->states_num / (int)pow(2, index);
    int new_index = 0;
    Ket new_state = ket_calloc(qupair->states_num);
    if (!complex_is_zero(alpha)) {
        for (int times = 0; times < qupair->states_num / gap; times++) {
            for (int i = 0; i < gap; ++i) {
                int offset = times * 2 * gap;
                Complex old_value = ket_get(qupair->state, offset + i);
                ket_set(new_state, new_index, complex_div(old_value, alpha));
                ++new_index;
            }
        }
    }
    else {
        for (int times = 0; times < qupair->states_num / gap; times++) {
            for (int i = 0; i < gap; ++i) {
                int offset = times * 2 * gap;
                Complex old_value = ket_get(qupair->state, offset + gap + i);
                ket_set(new_state, new_index, complex_div(old_value, beta));
                ++new_index;
            }
        }
    }
    ket_free(qupair->state);
    qupair->state = new_state;

    check_qupair(qupair);

    return 0;
}

static int attach_qubit_to_qureg(Qubit* qubit, Qureg* qureg, Complex alpha, Complex beta) {
    Qupair** new_qupairs = (Qupair**)pointer_memory_get(qureg->qupairs_num + 1);
    for (int i = 0; i < qureg->qupairs_num; ++i) {
        new_qupairs[i] = qureg->qupairs[i];
    }
    pointer_memory_return(qureg->qupairs, qureg->qupairs_num);
    ++qureg->qupairs_num;
    qureg->qupairs = new_qupairs;
    qureg->qupairs[qureg->qupairs_num - 1] = qupair_memory_get(1);
    initialize_qupair_with_single_qubit(qureg->qupairs[qureg->qupairs_num - 1], qureg, qureg->qupairs_num - 1, qubit, alpha, beta);

    return 0;
}

static int disentangle(Qupair* qupair, int separable_index, int is_normal, Complex ratio) {
    // calculate the probability amplitude according to the ratio
    Complex alpha, beta;
    calculate_probamp(ratio, is_normal, &alpha, &beta);

    // Update the newly-disentangled qubit
    Qubit* separable_qubit = qupair->qureg->qubits[qupair->qubits_indices[separable_index]];
    separable_qubit->entangled = 0;
    qubit_set_probamp(separable_qubit, alpha, beta);

    // Update the qupair
    detach_qubit_from_qupair(qupair, separable_index, alpha, beta);

    // Create a new qupair with single qubit, and attach it to the qureg
    attach_qubit_to_qureg(separable_qubit, separable_qubit->qureg, alpha, beta);

    return 0;
}

int check_entanglement(Qupair* qupair) {
    while (1) {
        // Qupair with single qubit means it's not entangled
        if (qupair->states_num == 2) {
            Qubit* left_qubit = qupair->qureg->qubits[qupair->qubits_indices[0]];
            left_qubit->entangled = 0;
            ket_set(left_qubit->state, 0, ket_get(qupair->state, 0));
            ket_set(left_qubit->state, 1, ket_get(qupair->state, 1));
            break;
        }

        int is_normal;
        Complex ratio;
        int separable_index = _check_entanglement(qupair, &ratio, &is_normal);

        if (separable_index != -1) {
            disentangle(qupair, separable_index, is_normal, ratio);
            // There maybe other unentangled qubits, so go on to check again
        }
        else {
            // All qubits in this qupair are entangled
            break;
        }
    }

    return 0;
}

int ket_positions_swap(Ket v, int sigs_num, int sig1, int sig2) {
    assert(sig1 < sigs_num && sig2 < sigs_num);
    assert(pow(2, sigs_num) == v.size);
    if (sig1 == sig2) {
        return 0;
    }
    int* swapped = int_memory_get(v.size);
    for (int i = 0; i < v.size; ++i) {
        swapped[i] = 0;
    }
    for (int i = 0; i < v.size; ++i) {
        // sigs are high-to-low
        if (swapped[i]) {
            continue;
        }
        int sig1_pos = (i >> (sigs_num - 1 - sig1)) & 0x1;
        int sig2_pos = (i >> (sigs_num - 1 - sig2)) & 0x1;
        if (sig1_pos == sig2_pos) {
            // don't need to swap
            swapped[i] = 1;
            continue;
        }
        int target_i = 0;
        for (int j = 0; j < sigs_num; ++j) {
            if (j == sig1) {
                target_i |= sig2_pos << (sigs_num - 1 - j);
            }
            else if (j == sig2) {
                target_i |= sig1_pos << (sigs_num - 1 - j);
            }
            else {
                target_i |= i & (0x1 << (sigs_num - 1 - j));
            }
        }
        Complex temp = ket_get(v, i);
        ket_set(v, i, ket_get(v, target_i));
        ket_set(v, target_i, temp);
        swapped[i] = 1;
        swapped[target_i] = 1;
    }
    int_memory_return(swapped, v.size);

    return 0;
}

Ket combine_qubits(Qubit** qubits, int qubits_num) {
    Ket combined_inputs;
    if (qubits_num == 1) {
        combined_inputs = ket_from_ket(qubits[0]->state);
    }
    else {
        combined_inputs = ket_Kronecker_product(qubits[0]->state, qubits[1]->state);
        for (int i = 2; i < qubits_num; ++i) {
            Ket temp = ket_Kronecker_product(combined_inputs, qubits[i]->state);
            ket_free(combined_inputs);
            combined_inputs = temp;
        }
    }

    return combined_inputs;
}

Ket combine_qupairs(Qupair** qupairs, int qupairs_num) {
    Ket combined_inputs;
    if (qupairs_num == 1) {
        combined_inputs = ket_from_ket(qupairs[0]->state);
    }
    else {
        combined_inputs = ket_Kronecker_product(qupairs[0]->state, qupairs[1]->state);
        for (int i = 2; i < qupairs_num; ++i) {
            Ket temp = ket_Kronecker_product(combined_inputs, qupairs[i]->state);
            ket_free(combined_inputs);
            combined_inputs = temp;
        }
    }

    return combined_inputs;
}

int get_exclusive_qupairs_from_qubits(Qubit** qubits, int qubits_num, Qupair** qupairs) {
    int qupairs_num = 0;
    for (int i = 0; i < qubits_num; ++i) {
        if (qubits[i]->entangled) {
            // Better use a 'set'
            int exist = 0;
            for (int j = 0; j < qupairs_num; ++j) {
                if (qupairs[j] == qubits[i]->qupair) {
                    exist = 1;
                    break;
                }
            }
            if (!exist) {
                qupairs[qupairs_num] = qubits[i]->qupair;
                ++qupairs_num;
            }
        }
    }

    return qupairs_num;
}

int get_entangled_qubits_num(Qupair** qupairs, int qupairs_num) {
    int num = 0;
    for (int i = 0; i < qupairs_num; ++i) {
        num += qupairs[i]->qubits_num;
    }

    return num;
}

int get_separable_qubits_num(Qubit** qubits, int qubits_num) {
    int num = 0;
    for (int i = 0; i < qubits_num; ++i) {
        if (!qubits[i]->entangled) {
            ++num;
        }
    }

    return num;
}

int unfold_qupairs_and_qubits(int* index_array, Qupair** qupairs, int qupairs_num, Qubit** qubits, int qubits_num) {
    int index = 0;
    for (int i = 0; i < qupairs_num; ++i) {
        for (int j = 0; j < qupairs[i]->qubits_num; ++j) {
            index_array[index] = qupairs[i]->qubits_indices[j];
            ++index;
        }
    }
    for (int i = 0; i < qubits_num; ++i) {
        if (!qubits[i]->entangled) {
            index_array[index] = qubits[i]->index;
            ++index;
        }
    }

    return 0;
}

int get_separable_qubits(Qubit** separable_qubits, Qubit** qubits, int qubits_num) {
    int separable_qubits_index = 0;
    for (int i = 0; i < qubits_num; ++i) {
        if (!qubits[i]->entangled) {
            separable_qubits[separable_qubits_index] = qubits[i];
            ++separable_qubits_index;
        }
    }

    return 0;
}
