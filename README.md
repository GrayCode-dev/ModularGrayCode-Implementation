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

The MPI implementation follows the same partitioning strategy as the OpenMP version, distributing codewords among processes. Each process independently generates its assigned subset without requiring inter-process communication and can store the results in a separate file. If the total number of generated codewords is needed, each process returns its local count, which can then be collected using MPI_Gather.

### Parallelization of recursive algorithm

The considered recursive algorithm generates the $m$-ary modular Gray code using recursion to emulate nested loops. Similar to the subset-generation algorithm, modulo operations are optimized with a simple **if statement**, reducing computational overhead. The parameter $k$ determines the last position to be modified during the generation process, while recursive calls control the coordinate changes. The algorithm starts from the vector $(0…01)$ with an initial value $k=1$.

```cpp
void Gray_rec(int g[],int k, int n, int m, int &br){
   if (k <= n) {
   	Gray_rec(g,k + 1, n, m,br);
   	for (int i = 1; i <= m - 1; i++) {
   	  if (g[k] == (m - 1)) g[k] = 0;
   	  else  (g[k]++);
   	  br++;
   	  Gray_rec(g,k + 1, n, m,br);
} } }
```
This parallelization approach fixes the first $nfix$ coordinates of the codewords and treats the generation of the remaining $n−nfix$ coordinates as independent computational tasks. Since there are $m^{nfix}$ possible combinations of fixed coordinates, the algorithm creates the same number of independent tasks that can be executed concurrently. In the OpenMP implementation, each task receives its own copy of the fixed coordinates, generates the corresponding subcode independently, and updates the total number of generated codewords through a **critical** section when finished. This task-based strategy enables efficient parallel execution on shared-memory systems.

```cpp
struct stg{int g[32];}sg;
int taskID = 0, total = 0;
void GrayTask(stg& sg, int k, int& nfix,
int& n, int& m) {
   if (k <= nfix) {
   	GrayTask(sg, k + 1, nfix, n, m);
   	for (int i = 1; i <= m - 1; i++) {
   	  if (sg.g[k] == (m - 1)) sg.g[k] = 0;
   	  else  (sg.g[k]++);
   	  taskID++;
# pragma omp task firstprivate(sg,taskID)
{
      sg.g[nfix + 1] = ((m - 1)*taskID)%m;
      int br = 0;
      Gray_rec(sg.g, nfix + 1,n, m, br);
# pragma omp critical
{total = total + br + 1; }
}
      GrayTask(sg, k + 1, nfix, n, m);
      }
  }
}

void mainTask(int nfix, int n, int m) {
	for (int i = 0; i <= n; i++) {
	   sg.g[i] = 0;
	   }
	   total++;
#pragma omp parallel num_threads(4)
{
#pragma omp single nowait
{
      int br = 0;
      for (int i = 0; i <= n; i++) {
      	sg.g[i] = 0;
      	}
      Gray_rec(sg.g, nfix + 1, n,
         m, br);
# pragma omp critical
{total = total + ct; }
      sg.g[nfix + 1] = m - 1;
      GrayTask(sg, 1, nfix, n, m);
}}}
```

The MPI implementation uses a **Master–Worker** strategy, where the master process dynamically distributes independent computational tasks to available worker processes. Each task corresponds to a fixed set of $nfix$ coordinates and generates the remaining codewords independently. Workers request new tasks upon completion, enabling efficient load balancing. The number of generated codewords can be collected from all processes through MPI collective communication.

### Experimental results

We present preliminary experimental results used to evaluate the efficiency of the proposed parallel implementations. The experimental results give the execution times in seconds. All computations are executed on Ubuntu 18.04 operation system with 2 CPUs Intel Xeon E5-2620v4 each with 8 cores and 16 threads. The MPI implementation is compiled with mpicc wrapper compiler of GNU Compiler Collection (GCC) version 7.5.0, while the OpenMP implementation uses the standard gcc compiler (version 7.5.0). The presented execution times are for implementations that calculate the total number of codewords without saving the words themselves.

Table 1 presents the execution times of the optimized sequential algorithms. The results show that the recursive algorithm performs better for smaller code lengths, while its efficiency decreases as $n$ increases due to the growing number of recursive calls. For larger lengths, the subset-generation algorithm becomes more efficient; for example, at $n=15$ the recursive version is approximately 41% slower. For smaller lengths, the performance difference between the two algorithms is below 10%.

| Algorithm  |m=5, n=15|m=11, n=10|m=29, n=7|m=31, n=7| 
|------------|---------|----------|---------|----------|
| Recursive  |  96,64  |  69,57   |  46,58  |  73,95   | 
| Subsets    |  68,17  |  75,03   |  47,71  |  76,29   | 
**Table 1**. Execution times using sequential implementations with optimized algorithms(sec.)

Tables 2 and 3 present the execution times of the parallel implementations for m=5,n=15 and m=11,n=10, respectively. In both cases, increasing the number of threads/processes reduces the execution time, although the speedup becomes less pronounced at higher levels of parallelism. The MPI implementations consistently outperform the OpenMP versions, with the MPI(GC) approach achieving the best execution times across all tested configurations. For OpenMP, the task-based implementation generally performs better than or similarly to the parallel-for version, particularly for larger workloads. Overall, the results demonstrate that MPI provides the highest scalability and efficiency, while both MPI and OpenMP benefit from increased parallelism.

| Algorithm    |  P=4  |  P=8  |  P=16  |  P=32  | 
|--------------|-------|-------|--------|--------|
| MPI(MW)      | 37,47 | 18,30 |  10,20 |  9,554 | 
| MPI(GC)      | 18,29 | 9,640 |  6,133 |  4,597 | 
| OpenMP(task) | 32,27 | 19,54 |  17,76 |  13,85 | 
| OpenMP(PF)   | 32,07 | 17,35 |  14,10 |  8,789 | 

**Table 2**. Execution times for $m=5$ $n=15$ $nfix=2$

| Algorithm    |  P=4  |  P=8  |  P=16  |  P=32  | 
|--------------|-------|-------|--------|--------|
| MPI(MW)      | 24,43 | 11,54 |  6,435 |  5,540 | 
| MPI(GC)      | 17,30 | 9,171 |  5,635 |  4,033 | 
| OpenMP(task) | 23,43 | 12,20 |  9,863 |  6,774 | 
| OpenMP(PF)   | 24,27 | 13,05 |  10,71 |  6,420 | 

**Table 3**. Execution times for $m=11$ $n=10$ $nfix=2$}

##### Related Publication

This repository accompanies research on the parallel generation of modular m-ary Gray codes. It contains the source code of the MPI and OpenMP implementations described in the corresponding publication.

More details about the implemented algorithms, parallelization strategies, performance analysis, and experimental results can be found in the following paper:

**D. Bikov, M. Pashinska-Gadzheva, and I. Bouyukliev**, *About Parallel Implementation of m-ary Gray Code*.
