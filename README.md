# Parallel implementation of m-ary Gray Code

This repository contains source code for parallel implementations of algorithms for generating modular m-ary Gray codes.

Two different parallelization approaches are included. The first approach uses an unranking function to partition the set of all codewords into independent subsets that can be generated concurrently. The second approach is based on a direct parallelization of a recursive generation algorithm.

The repository provides implementations based on both MPI (Message Passing Interface) and OpenMP. These implementations enable the evaluation and comparison of different parallelization strategies for recursive and non-recursive algorithms.

The source code is written in C/C++ and is intended for experimentation, performance evaluation, and further research in parallel combinatorial generation.


## Introduction

An ordering of all vectors in the $n$-dimensional vector space over $F_2$ is called a \emph{binary Gray code} if the Hamming distance between any two consecutive codewords is exactly one.

We can define $m$-ary Gray code as ordering of the vectors in $\Z_m^n$, where each pair of adjacent vector differs in one position. Two of the more popular $m$-ary Gray codes are the \emph{reflected} and the \emph{modular} Gray codes.



### Related Publication

This repository accompanies research on the parallel generation of modular m-ary Gray codes. It contains the source code of the MPI and OpenMP implementations described in the corresponding publication.

More details about the implemented algorithms, parallelization strategies, performance analysis, and experimental results can be found in the following paper:

**D. Bikov, M. Pashinska-Gadzheva, and I. Bouyukliev**, *About Parallel Implementation of m-ary Gray Code*.
