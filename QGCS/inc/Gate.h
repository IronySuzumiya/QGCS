#ifndef _Gate_H
#define _Gate_H

#include "Qubit.h"

extern Matrix H;

void apply_H(Qubit* qubit);
void apply_H_dagger(Qubit* qubit);
void apply_X(Qubit* qubit);
void apply_X_dagger(Qubit* qubit);
void apply_R(Qubit* qubit, float phi);
void apply_R_dagger(Qubit* qubit, float phi);
void apply_Z(Qubit* qubit);
void apply_Z_dagger(Qubit* qubit);
void apply_CNOT_on_int(Qubit** qubits, int qubits_num, int controlled);
void apply_CNOT_dagger_on_int(Qubit** qubits, int qubits_num, int controlled);
void apply_CNOT(Qubit** qubits, int qubits_num);
void apply_CNOT_dagger(Qubit** qubits, int qubits_num);
void apply_CR_on_int(Qubit** qubits, int qubits_num, int controlled, float phi);
void apply_CR_dagger_on_int(Qubit** qubits, int qubits_num, int controlled, float phi);
void apply_CR(Qubit** qubits, int qubits_num, float phi);
void apply_CR_dagger(Qubit** qubits, int qubits_num, float phi);
void apply_CZ_on_int(Qubit** qubits, int qubits_num, int controlled);
void apply_CZ_dagger_on_int(Qubit** qubits, int qubits_num, int controlled);
void apply_CZ(Qubit** qubits, int qubits_num);
void apply_CZ_dagger(Qubit** qubits, int qubits_num);
void apply_PauliZ_M(Qubit* qubit);

void apply_P(Qubit** qubits, int qubits_num);
void apply_P_dagger(Qubit** qubits, int qubits_num);

#endif
