#include <stdlib.h>
#include <string.h>

#include "lfa.h"

/* An implementation of the Knuth shuffling algorithm.  This can be done
   iteratively for paranoia purposes if the random numbers aren't very
   random.  Shuffles list of n elements of size s. */
void shuffle (struct rng_state *rs,
              void *list, void *tmp,
              size_t n, size_t s) {
  size_t felem, relem, nleft;
  unsigned char *array = list;

  /* For each element in the pack... */
  for(felem = 0, nleft = n; nleft > 0; felem++, nleft--) {
    /* Choose an element from the part of the pack we have not yet
       passed through. */
    relem = felem + rng_fetch_one_uniform(rs) * nleft;
    
    /* Swap them */
    memcpy(tmp, array + felem * s, s);
    memcpy(array + felem * s, array + relem * s, s);
    memcpy(array + relem * s, tmp, s);
  }
}

