#include <stdio.h>

void print_arr(int *a, int size) {
  while (size) {
    printf("%d ", a[size-1]);
    size--;
  }
  printf("\n");
}

int main() {
  int N = 10;
  int a[10];
  for (int i=0; i < N; i++)
    a[i] = i;
  print_arr(a, N);
  return 0;
}

