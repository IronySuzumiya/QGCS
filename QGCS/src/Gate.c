#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <memory.h>

#include "Gate.h"
#include "Tensor.h"

static gsl_vector_complex* matrix_vector_complex_mul(gsl_matrix_complex* m, gsl_vector_complex* v) {
    assert(m->size2 == v->size);
    gsl_vector_complex* temp = gsl_vector_complex_calloc(v->size);
    gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(1.0, 0), m, v, gsl_complex_rect(0, 0), temp);
    return temp;
}

static gsl_matrix_complex* matrix_matrix_complex_mul(gsl_matrix_complex* m1, gsl_matrix_complex* m2) {
    assert(m1->size2 == m2->size1);
    gsl_matrix_complex* temp = gsl_matrix_complex_calloc(m1->size1, m2->size2);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1.0, 0), m1, m2, gsl_complex_rect(0, 0), temp);
    return temp;
}

static bool complex_is_zero(gsl_complex v) {
    return GSL_REAL(v) > -1e-5 && GSL_REAL(v) < 1e-5 && GSL_IMAG(v) > -1e-5 && GSL_IMAG(v) < 1e-5;
}

static bool complex_equal(gsl_complex a, gsl_complex b) {
    return complex_is_zero(gsl_complex_sub(a, b));
}

static int _check_entanglement(Qupair* qupair) {
    struct {
        enum { unassigned, normal, infinite } type;
        gsl_complex value;
    } rate = { 0 };
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
        for (int times = 0; times < states_num / 2 / gap; times++) {
            bool entangled = false;
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
                        rate.type = infinite;
                    }
                    else {
                        rate.type = normal;
                        rate.value = gsl_complex_div(a, b);
                    }
                }
                else {
                    if (a_zero && b_zero) {
                        //always equal
                    }
                    else if (b_zero) {
                        if (rate.type != infinite) {
                            entangled = true;
                            break;
                        }
                    }
                    else {
                        if (rate.type != normal || !complex_equal(rate.value, gsl_complex_div(a, b))) {
                            entangled = true;
                            break;
                        }
                    }
                }
            }
            if (!entangled) {
                unentangled_index = index;
                goto unentangling;
            }
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
        old_qureg->qupairs[old_qureg->qupairs_num - 1] = (Qupair*)malloc(sizeof(Qupair));
        Qupair* new_qupair = old_qureg->qupairs[old_qureg->qupairs_num - 1];
        new_qupair->index = old_qureg->qupairs_num - 1;
        new_qupair->qureg = old_qureg;
        new_qupair->qubits_num = 1;
        new_qupair->qubits_indices = (int*)malloc(sizeof(int));
        new_qupair->qubits_indices[0] = unentangled_qubit->index;
        new_qupair->states_num = 2;
        new_qupair->state = gsl_vector_complex_calloc(2);
        gsl_vector_complex_set(new_qupair->state, 0, gsl_complex_rect(1.0, 0));
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

static void check_entanglement(Qupair* qupair) {
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
            }
            for (int j = 0; j < elements_per_cell / 2; ++j) {
                new_order[new_index] = (2 * i) * elements_per_cell / 2 + j;
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

static void apply_CNOT(Qubit** qubits, int qubits_num) {
    // For Simplicity
    Qureg* qureg = qubits[0]->qureg;
    for (int i = 1; i < qubits_num; ++i) {
        assert(qubits[i]->qureg == qureg);
    }

    // As usual, check entanglement first
    bool all_product = true;
    for (int i = 0; i < qubits_num; ++i) {
        if (qubits[i]->entangled) {
            all_product = false;
            break;
        }
    }

    // If all qubits are not entangled
    if (all_product) {
        gsl_vector_complex* combined_inputs = Kronecker_product_vv(qubits[0]->state, qubits[1]->state);
        for (int i = 2; i < qubits_num; ++i) {
            gsl_vector_complex* temp = Kronecker_product_vv(combined_inputs, qubits[i]->state);
            gsl_vector_complex_free(combined_inputs);
            combined_inputs = temp;
        }
        int states_num = combined_inputs->size;
        gsl_complex temp = gsl_vector_complex_get(combined_inputs, states_num - 2);
        gsl_vector_complex_set(combined_inputs, states_num - 2, gsl_vector_complex_get(combined_inputs, states_num - 1));
        gsl_vector_complex_set(combined_inputs, states_num - 1, temp);

        // Delete thoses old qupairs that these qubits belong to, and then make them into a new qupair
        int new_qupairs_num = qureg->qupairs_num - qubits_num + 1;
        Qupair** new_qupairs = (Qupair**)malloc(sizeof(Qupair*) * new_qupairs_num);
        // TODO
    }
    else {
        // Find all entangled qubits in the qureg
        bool* entangled = (int*)malloc(sizeof(bool) * qureg->qubits_num);
        memset(entangled, 0, sizeof(bool) * qureg->qubits_num);
        for (int i = 0; i < qubits_num; ++i) {
            if (qubits[i]->entangled) {
                for (int j = 0; j < qubits[i]->qupair->qubits_num; ++j) {
                    entangled[qubits[i]->qupair->qubits_indices[j]] = true;
                }
            }
        }
        // Check if all these entangled qubits are in the inputs of this gate
        for (int i = 0; i < qureg->qubits_num; ++i) {
            if (entangled[i]) {
                bool not_in = true;
                for (int i = 0; i < qubits_num; ++i) {
                    if (qubits[i]->index == i) {
                        not_in = false;
                        break;
                    }
                }
                assert(!not_in, "Just for now. The code will be over-complicated in this situation.");
            }
        }
        //   Here maybe some shuffles applied to qubits, because some of them are in their own qupairs
        // in which there are some spaces between their ids or they have different order compared to
        // that in the inputs.
        //   Quite tricky to deal with.
        // TODO
    }
}

static void apply_CNOT_dagger(Qubit** qubits, int qubits_num) {
    apply_CNOT(qubits, qubits_num);
}

ControlledGate CNOT = {
    apply_CNOT,
    apply_CNOT_dagger
};

int gate_init() {
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

    return 0;
}
