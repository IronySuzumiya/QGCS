# QGCS (Quantum Graph Computing Simulator)

A highly-optimized quantum simulator written in C, aimed at accelerating graph algorithms, planned to be deployed on a FPGA using HLS in the future.

A quantum minimum searching algorithm has been implemented, based on Grover's algorithm, serving as a component in SSSP.

To fit in OpenCL's specification, the dependency on GSL has been removed, and a simple memory allocation model also has been implemented.

The rest of code transplantation is still being worked on.
