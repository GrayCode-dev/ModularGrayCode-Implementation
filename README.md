# Parallel implementation of m-ary Gray Code

This repository contains source code for parallel implementations of algorithms for generating modular m-ary Gray codes.

Two different parallelization approaches are included. The first approach uses an unranking function to partition the set of all codewords into independent subsets that can be generated concurrently. The second approach is based on a direct parallelization of a recursive generation algorithm.

The repository provides implementations based on both MPI (Message Passing Interface) and OpenMP. These implementations enable the evaluation and comparison of different parallelization strategies for recursive and non-recursive algorithms.

The source code is written in C/C++ and is intended for experimentation, performance evaluation, and further research in parallel combinatorial generation.





### Related Publication

This repository accompanies the research work on parallel generation of modular m-ary Gray codes. The paper presents the underlying algorithms, parallelization strategies, experimental results, and a comparison of MPI and OpenMP implementations.
