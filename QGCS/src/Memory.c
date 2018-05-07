#include <assert.h>
#include <stdio.h>

#include "Memory.h"

typedef struct _complex_memory_block {
    Complex* head;
    Complex* tail;
    int size;
    struct _complex_memory_block* last;
    struct _complex_memory_block* next;
} ComplexMemoryBlock;

typedef struct _int_memory_block {
    int* head;
    int* tail;
    int size;
    struct _int_memory_block* last;
    struct _int_memory_block* next;
} IntMemoryBlock;

typedef struct _pointer_memory_block {
    void** head;
    void** tail;
    int size;
    struct _pointer_memory_block* last;
    struct _pointer_memory_block* next;
} PointerMemoryBlock;

typedef struct _qubit_memory_block {
    Qubit* head;
    Qubit* tail;
    int size;
    struct _qubit_memory_block* last;
    struct _qubit_memory_block* next;
} QubitMemoryBlock;

typedef struct _qupair_memory_block {
    Qupair* head;
    Qupair* tail;
    int size;
    struct _qupair_memory_block* last;
    struct _qupair_memory_block* next;
} QupairMemoryBlock;

typedef struct _qureg_memory_block {
    Qureg* head;
    Qureg* tail;
    int size;
    struct _qureg_memory_block* last;
    struct _qureg_memory_block* next;
} QuregMemoryBlock;

static ComplexMemoryBlock complex_memory_block[10000];
static IntMemoryBlock int_memory_block[1000];
static PointerMemoryBlock pointer_memory_block[500];
static QubitMemoryBlock qubit_memory_block[50];
static QupairMemoryBlock qupair_memory_block[50];
static QuregMemoryBlock qureg_memory_block[5];

static ComplexMemoryBlock* complex_memory_block_head = complex_memory_block;
static ComplexMemoryBlock* complex_memory_block_tail = &complex_memory_block[1];
static IntMemoryBlock* int_memory_block_head = int_memory_block;
static IntMemoryBlock* int_memory_block_tail = &int_memory_block[1];
static PointerMemoryBlock* pointer_memory_block_head = pointer_memory_block;
static PointerMemoryBlock* pointer_memory_block_tail = &pointer_memory_block[1];
static QubitMemoryBlock* qubit_memory_block_head = qubit_memory_block;
static QubitMemoryBlock* qubit_memory_block_tail = &qubit_memory_block[1];
static QupairMemoryBlock* qupair_memory_block_head = qupair_memory_block;
static QupairMemoryBlock* qupair_memory_block_tail = &qupair_memory_block[1];
static QuregMemoryBlock* qureg_memory_block_head = qureg_memory_block;
static QuregMemoryBlock* qureg_memory_block_tail = &qureg_memory_block[1];

static Complex complex_memory[10000000];
static int int_memory[5000];
static void* pointer_memory[500];
static Qubit qubit_memory[100];
static Qupair qupair_memory[100];
static Qureg qureg_memory[10];

#define _memory_init(blocks, num_block, memory, memory_size) \
    do {                                            \
        blocks[0].head = memory;                    \
        blocks[0].tail = &memory[memory_size];      \
        blocks[0].size = memory_size;               \
        for (int i = 0; i < num_block - 1; ++i) {   \
            blocks[i].next = &blocks[i + 1];        \
            blocks[i + 1].last = &blocks[i];        \
        }                                           \
        blocks[0].last = 0;                         \
        blocks[num_block - 1].next = 0;             \
    } while(0)

#define _block_delete(block_to_delete, block_tail) \
    do {                                                            \
        block_to_delete->head = 0;                                  \
        block_to_delete->tail = 0;                                  \
        block_to_delete->size = 0;                                  \
        if (block_to_delete->last) {                                \
            block_to_delete->last->next = block_to_delete->next;    \
        }                                                           \
        if (block_to_delete->next) {                                \
            block_to_delete->next->last = block_to_delete->last;    \
        }                                                           \
        block_tail->last->next = block_to_delete;                   \
        block_to_delete->last = block_tail->last;                   \
        block_to_delete->next = block_tail;                         \
        block_tail->last = block_to_delete;                         \
        block_tail = block_to_delete;                               \
    } while(0)

#define _block_merge_head(block_head, block_tail, block_to_merge, ptr) \
    do {                                                    \
        ptr = block_head;                                   \
        while (ptr != block_tail) {                         \
            if (ptr->tail == block_to_merge->head) {        \
                ptr->tail = block_to_merge->tail;           \
                ptr->size += block_to_merge->size;          \
                if (block_to_merge == block_head) {         \
                    block_head = block_to_merge->next;      \
                }                                           \
                _block_delete(block_to_merge, block_tail);  \
                break;                                      \
            }                                               \
            ptr = ptr->next;                                \
        }                                                   \
    } while(0)

#define _block_merge_tail(block_head, block_tail, block_to_merge, ptr) \
    do {                                                    \
        ptr = block_head;                                   \
        while (ptr != block_tail) {                         \
            if (ptr->head == block_to_merge->tail) {        \
                ptr->head = block_to_merge->head;           \
                ptr->size += block_to_merge->size;          \
                if (block_to_merge == block_head) {         \
                    block_head = block_to_merge->next;      \
                }                                           \
                _block_delete(block_to_merge, block_tail);  \
                break;                                      \
            }                                               \
            ptr = ptr->next;                                \
        }                                                   \
    } while(0)

#define _memory_return(block_head, block_tail, memory, memory_size, ptr1, ptr2) \
    do {                                                                \
        assert(block_tail != 0);                                        \
        int mergable = 0;                                               \
        ptr1 = block_head;                                              \
        while (ptr1 != block_tail) {                                    \
            if (ptr1->head == &memory[memory_size]) {                   \
                ptr1->head = memory;                                    \
                ptr1->size += memory_size;                              \
                _block_merge_head(block_head, block_tail, ptr1, ptr2);  \
                mergable = 1;                                           \
                break;                                                  \
            }                                                           \
            else if (ptr1->tail == memory) {                            \
                ptr1->tail = &memory[memory_size];                      \
                ptr1->size += memory_size;                              \
                _block_merge_tail(block_head, block_tail, ptr1, ptr2);  \
                mergable = 1;                                           \
                break;                                                  \
            }                                                           \
            ptr1 = ptr1->next;                                          \
        }                                                               \
        if (!mergable) {                                                \
            block_tail->head = memory;                                  \
            block_tail->tail = &memory[memory_size];                    \
            block_tail->size = memory_size;                             \
            block_tail = block_tail->next;                              \
        }                                                               \
    } while(0)

#define _memory_get(block_head, block_tail, memory, memory_size, ptr) \
    do {                                                \
        ptr = block_head;                               \
        while (ptr != block_tail) {                     \
            if (ptr->size >= memory_size) {             \
                memory = ptr->head;                     \
                ptr->head = &ptr->head[memory_size];    \
                ptr->size -= memory_size;               \
                break;                                  \
            }                                           \
            ptr = ptr->next;                            \
        }                                               \
        assert(ptr != block_tail);                      \
        if (ptr->size == 0) {                           \
            if (ptr == block_head) {                    \
                block_head = ptr->next;                 \
            }                                           \
            _block_delete(ptr, block_tail);             \
        }                                               \
    } while(0)

void memory_init() {
    _memory_init(complex_memory_block, 10000, complex_memory, 10000000);
    _memory_init(int_memory_block, 1000, int_memory, 5000);
    _memory_init(pointer_memory_block, 500, pointer_memory, 500);
    _memory_init(qubit_memory_block, 50, qubit_memory, 100);
    _memory_init(qupair_memory_block, 50, qupair_memory, 100);
    _memory_init(qureg_memory_block, 5, qureg_memory, 10);
}

Complex* complex_memory_get(int size) {
    Complex* mem;
    ComplexMemoryBlock* ptr;
    /*if (size >= 4096) {
        print_memory_usage();
    }*/
    _memory_get(complex_memory_block_head, complex_memory_block_tail, mem, size, ptr);
    return mem;
}

void complex_memory_return(Complex* addr, int size) {
    ComplexMemoryBlock* ptr1;
    ComplexMemoryBlock* ptr2;
    _memory_return(complex_memory_block_head, complex_memory_block_tail, addr, size, ptr1, ptr2);
}

int* int_memory_get(int size) {
    int* mem;
    IntMemoryBlock* ptr;
    _memory_get(int_memory_block_head, int_memory_block_tail, mem, size, ptr);
    return mem;
}

void int_memory_return(int* addr, int size) {
    IntMemoryBlock* ptr1;
    IntMemoryBlock* ptr2;
    _memory_return(int_memory_block_head, int_memory_block_tail, addr, size, ptr1, ptr2);
}

void** pointer_memory_get(int size) {
    void** mem;
    PointerMemoryBlock* ptr;
    _memory_get(pointer_memory_block_head, pointer_memory_block_tail, mem, size, ptr);
    return mem;
}

void pointer_memory_return(void** addr, int size) {
    PointerMemoryBlock* ptr1;
    PointerMemoryBlock* ptr2;
    _memory_return(pointer_memory_block_head, pointer_memory_block_tail, addr, size, ptr1, ptr2);
}

Qubit* qubit_memory_get(int size) {
    Qubit* mem;
    QubitMemoryBlock* ptr;
    _memory_get(qubit_memory_block_head, qubit_memory_block_tail, mem, size, ptr);
    return mem;
}

void qubit_memory_return(Qubit* addr, int size) {
    QubitMemoryBlock* ptr1;
    QubitMemoryBlock* ptr2;
    _memory_return(qubit_memory_block_head, qubit_memory_block_tail, addr, size, ptr1, ptr2);
}

Qupair* qupair_memory_get(int size) {
    Qupair* mem;
    QupairMemoryBlock* ptr;
    _memory_get(qupair_memory_block_head, qupair_memory_block_tail, mem, size, ptr);
    return mem;
}

void qupair_memory_return(Qupair* addr, int size) {
    QupairMemoryBlock* ptr1;
    QupairMemoryBlock* ptr2;
    _memory_return(qupair_memory_block_head, qupair_memory_block_tail, addr, size, ptr1, ptr2);
}

Qureg* qureg_memory_get(int size) {
    Qureg* mem;
    QuregMemoryBlock* ptr;
    _memory_get(qureg_memory_block_head, qureg_memory_block_tail, mem, size, ptr);
    return mem;
}

void qureg_memory_return(Qureg* addr, int size) {
    QuregMemoryBlock* ptr1;
    QuregMemoryBlock* ptr2;
    _memory_return(qureg_memory_block_head, qureg_memory_block_tail, addr, size, ptr1, ptr2);
}

#define _print_memory_usage(memory, memory_size, block_head, block_tail, ptr) \
    do {                                                                        \
        printf("total:\n");                                                     \
        printf("  %p - %p  %d\n", memory, &memory[memory_size], memory_size);   \
        printf("free space:\n");                                                \
        ptr = block_head;                                                       \
        while (ptr != block_tail) {                                             \
            printf("  %p - %p  %d\n", ptr->head, ptr->tail, ptr->size);         \
            ptr = ptr->next;                                                    \
        }                                                                       \
    } while(0)

void print_memory_usage() {
    ComplexMemoryBlock* pcom;
    IntMemoryBlock* pint;
    PointerMemoryBlock* pp;
    QubitMemoryBlock* pqb;
    QupairMemoryBlock* pqp;
    QuregMemoryBlock* pqr;

    printf("complex memory:\n");
    _print_memory_usage(complex_memory, 10000000, complex_memory_block_head, complex_memory_block_tail, pcom);
    printf("\n");

    printf("int memory:\n");
    _print_memory_usage(int_memory, 5000, int_memory_block_head, int_memory_block_tail, pint);
    printf("\n");

    printf("pointer memory:\n");
    _print_memory_usage(pointer_memory, 500, pointer_memory_block_head, pointer_memory_block_tail, pp);
    printf("\n");

    printf("qubit memory:\n");
    _print_memory_usage(qubit_memory, 100, qubit_memory_block_head, qubit_memory_block_tail, pqb);
    printf("\n");

    printf("qupair memory:\n");
    _print_memory_usage(qupair_memory, 100, qupair_memory_block_head, qupair_memory_block_tail, pqp);
    printf("\n");

    printf("qureg memory:\n");
    _print_memory_usage(qureg_memory, 10, qureg_memory_block_head, qureg_memory_block_tail, pqr);
    printf("\n");

    printf("\n");
}
