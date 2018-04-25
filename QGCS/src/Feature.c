#include <memory.h>
#include <assert.h>
#include <math.h>

#include "Feature.h"
#include "Util.h"
#include "Tensor.h"

static int _check_entanglement(Qupair* qupair, gsl_complex* ratio, bool* is_normal) {
    enum { unassigned, normal, infinite } type;
    gsl_complex temp_value = gsl_complex_rect(0, 0);
    int separable_index = -1;
    int states_num = qupair->states_num;

    for (int gap = states_num / 2, index = 0; gap >= 1; gap /= 2, ++index) {
        type = unassigned;
        *ratio = gsl_complex_rect(0.0, 0.0);
        bool entangled = false;
        for (int times = 0; times < states_num / 2 / gap; times++) {
            for (int i = 0; i < gap; ++i) {
                int offset = times * 2 * gap;
                gsl_complex a = gsl_vector_complex_get(qupair->state, offset + i);
                gsl_complex b = gsl_vector_complex_get(qupair->state, offset + gap + i);
                bool a_zero = complex_is_zero(a);
                bool b_zero = complex_is_zero(b);
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
                        *ratio = gsl_complex_div(a, b);
                    }
                }
                else {
                    if (a_zero && b_zero) {
                        //always equal
                    }
                    else if (b_zero) {
                        if (type != infinite || !complex_equal(a, temp_value)) {
                            entangled = true;
                            break;
                        }
                    }
                    else {
                        if ((!complex_is_zero(temp_value) && a_zero && !complex_equal(b, temp_value))
                            || type != normal || !complex_equal(*ratio, gsl_complex_div(a, b))) {
                            entangled = true;
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

static int calculate_probamp(gsl_complex ratio, bool is_normal, gsl_complex* alpha, gsl_complex* beta) {
    gsl_complex one = gsl_complex_rect(1.0, 0.0);
    gsl_complex zero = gsl_complex_rect(0.0, 0.0);
    if (is_normal) {
        gsl_complex r2_plus_1 = gsl_complex_add(gsl_complex_mul(ratio, ratio), one);
        gsl_complex sqrt_1_over_r2_plus_1 = gsl_complex_sqrt(gsl_complex_div(one, r2_plus_1));
        *alpha = gsl_complex_mul(ratio, sqrt_1_over_r2_plus_1);
        *beta = sqrt_1_over_r2_plus_1;
    }
    else {
        *alpha = one;
        *beta = zero;
    }

    return 0;
}

static int detach_qubit_from_qupair(Qupair* qupair, int index, gsl_complex alpha, gsl_complex beta) {
    --qupair->qubits_num;
    qupair->states_num /= 2;
    // Avoid to use realloc()
    int* new_qubits_indices = (int*)malloc(sizeof(int) * qupair->qubits_num);
    memcpy(new_qubits_indices, qupair->qubits_indices, sizeof(int) * index);
    memcpy(&new_qubits_indices[index],
        &qupair->qubits_indices[index + 1], sizeof(int) * (qupair->qubits_num - index));
    free(qupair->qubits_indices);
    qupair->qubits_indices = new_qubits_indices;

    int gap = qupair->states_num / (int)pow(2, index);
    int new_index = 0;
    gsl_vector_complex* new_state = gsl_vector_complex_calloc(qupair->states_num);
    for (int times = 0; times < qupair->states_num / gap; times++) {
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
    gsl_vector_complex_free(qupair->state);
    qupair->state = new_state;

    return 0;
}

static int attach_qubit_to_qureg(Qubit* qubit, Qureg* qureg, gsl_complex alpha, gsl_complex beta) {
    ++qureg->qupairs_num;
    Qupair** new_qupairs = (Qupair**)malloc(sizeof(Qupair*) * qureg->qupairs_num);
    for (int i = 0; i < qureg->qupairs_num - 1; ++i) {
        new_qupairs[i] = qureg->qupairs[i];
    }
    free(qureg->qupairs);
    qureg->qupairs = new_qupairs;
    qureg->qupairs[qureg->qupairs_num - 1] = (Qupair*)malloc(sizeof(Qupair));
    initialize_qupair_with_single_qubit(qureg->qupairs[qureg->qupairs_num - 1], qureg, qureg->qupairs_num - 1, qubit, alpha, beta);

    return 0;
}

static int disentangle(Qupair* qupair, int separable_index, bool is_normal, gsl_complex ratio) {
    // calculate the probability amplitude according to the ratio
    gsl_complex alpha, beta;
    calculate_probamp(ratio, is_normal, &alpha, &beta);

    // Update the newly-disentangled qubit
    Qubit* separable_qubit = qupair->qureg->qubits[qupair->qubits_indices[separable_index]];
    separable_qubit->entangled = false;
    qubit_set_probamp(separable_qubit, alpha, beta);

    // Update the qupair
    detach_qubit_from_qupair(qupair, separable_index, alpha, beta);

    // Create a new qupair with single qubit, and attach it to the qureg
    attach_qubit_to_qureg(separable_qubit, separable_qubit->qureg, alpha, beta);

    return 0;
}

int check_entanglement(Qupair* qupair) {
    while (true) {
        // Qupair with single qubit means it's not entangled
        if (qupair->states_num == 2) {
            Qubit* left_qubit = qupair->qureg->qubits[qupair->qubits_indices[0]];
            left_qubit->entangled = false;
            gsl_vector_complex_memcpy(left_qubit->state, qupair->state);
            break;
        }

        bool is_normal;
        gsl_complex ratio;
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

int vector_complex_positions_swap(gsl_vector_complex* v, int sigs_num, int sig1, int sig2) {
    assert(sig1 < sigs_num && sig2 < sigs_num);
    assert(pow(2, sigs_num) == v->size);
    if (sig1 == sig2) {
        return 0;
    }
    bool* swapped = (bool*)malloc(sizeof(bool) * v->size);
    memset(swapped, 0, sizeof(bool) * v->size);
    for (int i = 0; i < (int)v->size; ++i) {
        // sigs are high-to-low
        if (swapped[i]) {
            continue;
        }
        int sig1_pos = (i >> (sigs_num - 1 - sig1)) & 0x1;
        int sig2_pos = (i >> (sigs_num - 1 - sig2)) & 0x1;
        if (sig1_pos == sig2_pos) {
            // don't need to swap
            swapped[i] = true;
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
        gsl_complex temp = gsl_vector_complex_get(v, i);
        gsl_vector_complex_set(v, i, gsl_vector_complex_get(v, target_i));
        gsl_vector_complex_set(v, target_i, temp);
        swapped[i] = true;
        swapped[target_i] = true;
    }
    free(swapped);

    return 0;
}

gsl_vector_complex* combine_qubits(Qubit** qubits, int qubits_num) {
    gsl_vector_complex* combined_inputs;
    if (qubits_num == 1) {
        combined_inputs = gsl_vector_complex_calloc(2);
        gsl_vector_complex_memcpy(combined_inputs, qubits[0]->state);
    }
    else {
        combined_inputs = Kronecker_product_vv(qubits[0]->state, qubits[1]->state);
        for (int i = 2; i < qubits_num; ++i) {
            gsl_vector_complex* temp = Kronecker_product_vv(combined_inputs, qubits[i]->state);
            gsl_vector_complex_free(combined_inputs);
            combined_inputs = temp;
        }
    }

    return combined_inputs;
}

gsl_vector_complex* combine_qupairs(Qupair** qupairs, int qupairs_num) {
    gsl_vector_complex* combined_inputs;
    if (qupairs_num == 1) {
        combined_inputs = gsl_vector_complex_calloc(qupairs[0]->states_num);
        gsl_vector_complex_memcpy(combined_inputs, qupairs[0]->state);
    }
    else {
        combined_inputs = Kronecker_product_vv(qupairs[0]->state, qupairs[1]->state);
        for (int i = 2; i < qupairs_num; ++i) {
            gsl_vector_complex* temp = Kronecker_product_vv(combined_inputs, qupairs[i]->state);
            gsl_vector_complex_free(combined_inputs);
            combined_inputs = temp;
        }
    }

    return combined_inputs;
}

int get_exclusive_qupairs_from_qubits(Qubit** qubits, int qubits_num, Qupair** qupairs) {
    int qupairs_num = 0;
    for (int i = 0; i < qubits_num; ++i) {
        if (qubits[i]->entangled) {
            // Better to use a 'set'
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
