#include <cstdio>

const int Maxn = 30;
int g[Maxn];

int a, b;

// print function
void print(int g[], int n, int t) {
    for (int j = 1; j <= n; j++)
        printf("%d", g[j]);
    printf(" T = %d\n", t);
}

void Gray_rec(int g[], int k, int n, int m, int &br) {
    if (k <= n) {
        Gray_rec(g, k + 1, n, m, br);

        for (int i = 1; i <= m - 1; i++) {
            if (g[k] == (m - 1))
                g[k] = 0;
            else
                (g[k]++);

            br++;
            Gray_rec(g, k + 1, n, m, br);
        }
    }
    // print Gray codewords
    else {
        print(g, n, br);
    }
    //=======================
}


void Gray_modular_unrank(int g[], int m, int n, int num) {
    int a[Maxn]={0};

    for (int i = 0; i < n; i++) {
        a[n - i] = num % m;
        num = num / m;
    }

    g[1] = a[1];

    for (int i = 2; i <= n; i++) {
        g[i] = (m - a[i - 1] + a[i]) % m;
    }
}

void Gray_subset(int m, int n, int A, int B){
   int helper[Maxn];
   for (int i = 0; i <= n; i++) {
      helper[i] = 0; }
   Gray_modular_unrank(g, m, n, a);
   int i = a, j=a, k, t = n, br = 0;
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
            g[k] = 0; t = k;}
         else {g[k] = (g[k] + 1); t = k;}
      }
      br++;
      print(g, n, br); // print next codeword
   } while (k > 0);
}


int main() {
    // setting up parameters
    int n = 3;      // word length
    int m = 3;      // alphabet size
    int br = 0;     // codewords counter


    for (int i = 0; i < Maxn; i++)
        g[i] = 0;

    //Call Recursive function
    Gray_rec(g, 1, n, m, br);

    printf("\nTotal transitions (Recursive) = %d\n\n", br);
    //=======================

    // setting up parameters for Subsets implementation
    a = 0;  // starting rank
    b = br;  // ending rank

    Gray_modular_unrank(g, m, n, a);
    print(g, n, a); //print first codeword

    //Call Subsets function
    Gray_subset(m, n, a, b);

    return 0;
}
