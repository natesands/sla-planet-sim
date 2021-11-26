#include <stdio.h>

int N;
void doubleN() {
  N = N*2;
}
int main() {
  N = 3;
  printf("%d\n",N);
  doubleN();
  printf("%d\n",N);
  return 0;
}

