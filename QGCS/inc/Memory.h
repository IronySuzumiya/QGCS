#ifndef _MEMORY_H
#define _MEMORY_H

#include "Qubit.h"

void memory_init(void);
Complex* complex_memory_get(int size);
void complex_memory_return(Complex* addr, int size);
int* int_memory_get(int size);
void int_memory_return(int* addr, int size);
void** pointer_memory_get(int size);
void pointer_memory_return(void** addr, int size);
Qubit* qubit_memory_get(int size);
void qubit_memory_return(Qubit* addr, int size);
Qupair* qupair_memory_get(int size);
void qupair_memory_return(Qupair* addr, int size);
Qureg* qureg_memory_get(int size);
void qureg_memory_return(Qureg* addr, int size);
void print_memory_usage(void);

#endif
