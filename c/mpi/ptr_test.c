#include <stdio.h>
int main() {

  int arr[10];
  for (int i=0; i < 10; i++)
    arr[i] = i;

  int *ptr;
  ptr = &arr[5];
  printf("%d\n", *ptr);

  return 0;
}



