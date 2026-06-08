# Parallel implementation of m-ary Gray Code

This repository contains source code for parallel implementations of algorithms for generating modular m-ary Gray codes.

Two different parallelization approaches are included. The first approach uses an unranking function to partition the set of all codewords into independent subsets that can be generated concurrently. The second approach is based on a direct parallelization of a recursive generation algorithm.

The repository provides implementations based on both MPI (Message Passing Interface) and OpenMP. These implementations enable the evaluation and comparison of different parallelization strategies for recursive and non-recursive algorithms.

The source code is written in C/C++ and is intended for experimentation, performance evaluation, and further research in parallel combinatorial generation.


### Introduction

An ordering of all vectors in the $n$-dimensional vector space over $F_2$ is called a $binary$ $Gray$ $code$ if the Hamming distance between any two consecutive codewords is exactly one.

We can define $m$-ary Gray code as ordering of the vectors in $Z_m^n$, where each pair of adjacent vector differs in one position. Two of the more popular $m$-ary Gray codes are the $reflected$ and the $modular$ Gray codes.

A modular $m$-ary Gray code can be defined as the ordering of $Z_m^n$, where the difference between two pairs of vectors is defined in terms of $directional$ $distance$. We denote directional distance as $\sigma(b,a) = (b-a)\pmod{m}$, where $a,b \in \Z_m$. In terms of vectors $u,v\in \Z_m^n$ it is defined as follows
$$d_D(u,v) = \sum_{i=1}^n\sigma(b_i - a_i),$$
where $u=(a_1,a_2,\dots,a_n), v = (b_1,b_2,\dots, b_n)$.
In the modular $m$-ary Gray code, the directional difference between a vector and its preceding is one.

The main problem with development of a parallel implementation of algorithms is to appropriately divide the computational work in such way that each process or thread would execute approximately equal amount of autonomous calculations. This division of work is dependent on the problem itself and the hardware architecture. Here we consider parallel implementations for architectures with shared and distributed memory. We use Open Multi-Processing (OpenMP) and Message Passing Interface (MPI) application programming interfaces for C/C++ for these architectures, respectively.

### Parallelization through generating subsets of the code

...

### Parallelization of recursive algorithm

...

### Experimental results

...

##### Related Publication

This repository accompanies research on the parallel generation of modular m-ary Gray codes. It contains the source code of the MPI and OpenMP implementations described in the corresponding publication.

More details about the implemented algorithms, parallelization strategies, performance analysis, and experimental results can be found in the following paper:

**D. Bikov, M. Pashinska-Gadzheva, and I. Bouyukliev**, *About Parallel Implementation of m-ary Gray Code*.
