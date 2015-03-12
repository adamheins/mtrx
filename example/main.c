
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "vctr.h"
#include "mtrx.h"

int main() {

  srand(time(NULL));

  matrix_t *A = mtrx_rnd(3, 3, 10);
  matrix_t *B = mtrx_rnd(3, 3, 10);
  matrix_t *C = mtrx_mult(A, B);

  printf("\nWelcome to the matrix.\n\n");
  mtrx_print(A);
  printf("\nmultiplied by\n\n");
  mtrx_print(B);
  printf("\nequals\n\n");
  mtrx_print(C);
  printf("\n");

  mtrx_destroy(A);
  mtrx_destroy(B);
  mtrx_destroy(C);
}
