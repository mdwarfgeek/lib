#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stddef.h>
#include <inttypes.h>

#include "lfa.h"

#define N 1024

int main (int argc, char *argv[]) {
  struct rng_state ds;
  double a[N];
  int i;

  rng_init(&ds, 42);

  rng_fetch_uniform(&ds, a, N);

  for(i = 0; i < N; i++)
    printf("%d %.20f\n", i, a[i]);

  return(0);
}
