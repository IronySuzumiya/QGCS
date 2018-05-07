#ifndef _FEATURE_H
#define _FEATURE_H

#include "Qubit.h"

int check_entanglement(Qupair* qupair);
int ket_positions_swap(Ket v, int sigs_num, int sig1, int sig2);
Ket combine_qubits(Qubit** qubits, int qubits_num);
Ket combine_qupairs(Qupair** qupairs, int qupairs_num);
int get_exclusive_qupairs_from_qubits(Qubit** qubits, int qubits_num, Qupair** qupairs);
int get_entangled_qubits_num(Qupair** qupairs, int qupairs_num);
int get_separable_qubits_num(Qubit** qubits, int qubits_num);
int unfold_qupairs_and_qubits(int* index_array, Qupair** qupairs, int qupairs_num, Qubit** qubits, int qubits_num);
int get_separable_qubits(Qubit** separable_qubits, Qubit** qubits, int qubits_num);

#endif
