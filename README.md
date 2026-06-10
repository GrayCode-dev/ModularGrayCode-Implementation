# Parallel implementation of m-ary Gray Code

This repository contains source code for parallel implementations of algorithms for generating modular m-ary Gray codes.

Two different parallelization approaches are included. The first approach uses an unranking function to partition the set of all codewords into independent subsets that can be generated concurrently. The second approach is based on a direct parallelization of a recursive generation algorithm.

The repository provides implementations based on both MPI (Message Passing Interface) and OpenMP. These implementations enable the evaluation and comparison of different parallelization strategies for recursive and non-recursive algorithms.

The source code is written in C/C++ and is intended for experimentation, performance evaluation, and further research in parallel combinatorial generation.


### Introduction

An ordering of all vectors in the $n$-dimensional vector space over $F_2$ is called a $binary$ $Gray$ $code$ if the Hamming distance between any two consecutive codewords is exactly one.

We can define $m$-ary Gray code as ordering of the vectors in $Z_m^n$, where each pair of adjacent vector differs in one position. Two of the more popular $m$-ary Gray codes are the $reflected$ and the $modular$ Gray codes.

A modular $m$-ary Gray code can be defined as the ordering of $Z_m^n$, where the difference between two pairs of vectors is defined in terms of $directional$ $distance$. We denote directional distance as $\sigma(b,a) = (b-a)\pmod{m}$, where $a,b \in Z_m$. In terms of vectors $u,v\in Z_m^n$ it is defined as follows
$$d_D(u,v) = \sum_{i=1}^n\sigma(b_i - a_i),$$
where $u=(a_1,a_2,\dots,a_n), v = (b_1,b_2,\dots, b_n)$.
In the modular $m$-ary Gray code, the directional difference between a vector and its preceding is one.

The main problem with development of a parallel implementation of algorithms is to appropriately divide the computational work in such way that each process or thread would execute approximately equal amount of autonomous calculations. This division of work is dependent on the problem itself and the hardware architecture. Here we consider parallel implementations for architectures with shared and distributed memory. We use Open Multi-Processing (OpenMP) and Message Passing Interface (MPI) application programming interfaces for C/C++ for these architectures, respectively.

### Parallelization through generating subsets of the code

The considered algorithm efficiently generates a subset of the modular $m$-ary Gray code $G(m,n)$ by producing only the codewords whose ranks belong to a specified interval $[A,B]$. It starts by applying an unranking procedure to obtain the codeword of rank $A$, and then sequentially generates the remaining $B−A$ codewords.

To improve performance, the implementation avoids repeated modulo division operations that are normally required to determine the coordinate change between consecutive codewords. Instead, it uses an auxiliary helper array, initialized from the m-ary representation of the starting rank A, together with a simple if statement to emulate the modulo calculations. This optimization efficiently tracks the coordinate updates needed for generating subsequent codewords, resulting in a faster generation of Gray code subsets.

```cpp
const int Maxn = 30;
int g[Maxn];
void Gray_modular_unrank(int g[], int m, int n, int a);

void Gray_subset(int m, int n, int A, int B){
   int helper[Maxn];
   for (int i = 0; i <= n; i++) {
   	helper[i] = 0; }
   Gray_modular_unrank(g, m, n, a);
   int i = a, j=a, k, t = n,br = 0;
   while (j > 0) {
     helper[t] = j % m;
     j = j / m; t--;}
   t = 0;
   do {
   	 if (br == b - a)  return;
   	 k = n;
   	 i = i + 1;
   	 j = i;
   	 helper[k]++;
   	 while (m == helper[k]) {
   	   helper[k] = 0;
   	   k = k - 1;
   	   helper[k]++;}
   	 if (k > 0) {
   	   if (g[k] == (m - 1)){
   	   	g[k] = 0;t = k;}
   	   else {g[k] = (g[k] + 1); t = k;}
   	 }
   	 br++;
   } while (k > 0);}
```
The considered parallel algorithm generates the complete $m$-ary Gray code with $A=0$ and $B=m^n-1$ by partitioning the range of codeword ranks among multiple threads or processes. Each thread independently generates a contiguous subset of codewords using the sequential generation algorithm. In an OpenMP implementation, the workload is distributed statically so that each thread computes approximately $m^n/P$ codewords, where $P$ is the number of threads. This approach enables efficient parallel generation of Gray code subsets and scalable execution on shared-memory systems.

```cpp
void GrayOpenMP_subset(int m, int n) {
   int ALLcw = 1;
   for (int i = 1; i <= n; i++) {
   	ALLcw = ALLcw * m;}
   	long long int P = 4;
   	long long int cwPerThread = ALLcw / P;
   	omp_set_num_threads(P);
#pragma omp parallel for
    for(int i = 0; i<P; i++){
    int A = i * cwPerThread,
     B = ((i + 1) * words_per_process) - 1;
    if (i == (P - 1)) B = total - 1;
    Gray_subset(m, n, A, B);
}}
```
To count the total number of generated codewords, a shared global counter can be updated inside an OpenMP $critical$ $section$. This ensures that only one thread modifies the counter at a time, preventing race conditions.

### Parallelization of recursive algorithm

...

### Experimental results

We present preliminary experimental results used to evaluate the efficiency of the proposed parallel implementations. The experimental results give the execution times in seconds. All computations are executed on Ubuntu 18.04 operation system with 2 CPUs Intel Xeon E5-2620v4 each with 8 cores and 16 threads. The MPI implementation is compiled with mpicc wrapper compiler of GNU Compiler Collection (GCC) version 7.5.0, while the OpenMP implementation uses the standard gcc compiler (version 7.5.0). The presented execution times are for implementations that calculate the total number of codewords without saving the words themselves.

...

##### Related Publication

This repository accompanies research on the parallel generation of modular m-ary Gray codes. It contains the source code of the MPI and OpenMP implementations described in the corresponding publication.

More details about the implemented algorithms, parallelization strategies, performance analysis, and experimental results can be found in the following paper:

**D. Bikov, M. Pashinska-Gadzheva, and I. Bouyukliev**, *About Parallel Implementation of m-ary Gray Code*.
