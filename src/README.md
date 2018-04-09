C++ implementation of a non-negative matrix factorization (NMF) algorithm.

# Introduction

This repository contains an implementation of alternating non-negative least 
squares (ANLS) with block-principal pivoting (BPP) for NMF. A partial GPU 
implementation (for level-3 matrix computations) is included.

# Compilation

To compile the serial version, simply compile the main_serial.cpp file 
(or use ./compile_serial).

To compile the parallel version, simply use 'make'.