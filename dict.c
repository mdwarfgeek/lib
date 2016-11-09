#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "lfa.h"

/* 32-bit FNV-1a hash */
#define HASH_INIT   2166136261
#define HASH_PRIME  16777619
#define HASH_ACCUM(h, b) (h) ^= (uint32_t) (b); (h) *= HASH_PRIME

static inline uint32_t hash_mem (void *str, uint32_t n) {
  uint32_t h, i;
  uint8_t *p;
  
  h = HASH_INIT;
  for(i = 0, p = (uint8_t *) str; i < n; i++, p++) {
    HASH_ACCUM(h, *p);
  }

  return(h);
}

int dict_init (struct dict *d, uint32_t nalloc) {
  if(nalloc > 0) {
    d->list
      = (struct dict_entry **) malloc(nalloc * sizeof(struct dict_entry *));
    if(!d->list)
      return(-1);
  }
  else
    d->list = (void *) NULL;

  d->nalloc = nalloc;
  d->nused = 0;

  return(0);
}

void dict_free (struct dict *d) {
  struct dict_entry *de, *pde;
  uint32_t ient;

  for(ient = 0; ient < d->nalloc; ient++) {
    if(d->list[ient]) {
      de = d->list[ient];
      while(de) {
        pde = de;
        de = de->next;

        free((void *) pde->key);
        free((void *) pde);
      }
    }
  }

  if(d->list)
    free((void *) d->list);

  d->nalloc = 0;
  d->nused = 0;
}

int dict_delete (struct dict *d,
                 void *key, uint32_t keylen) {
  uint32_t hash, ient;
  struct dict_entry *de, *ode;
  int rv = 0;

  hash = hash_mem(key, keylen);

  ient = hash % d->nalloc;

  de = d->list[ient];
  ode = (struct dict_entry *) NULL;

  while(de) {
    if(de->hash == hash &&
       de->keylen == keylen &&
       !memcmp(de->key, key, keylen)) {
      if(ode)
        ode->next = de->next;
      else
        d->list[ient] = (struct dict_entry *) NULL;

      free((void *) de->key);
      free((void *) de);

      d->nused--;

      rv = 1;

      break;
    }

    ode = de;
    de = de->next;
  }

  return(rv);
}

int dict_fetch (struct dict *d,
                void *key, uint32_t keylen,
                void **value) {
  uint32_t hash, ient;
  struct dict_entry *de;
  int rv = 0;

  hash = hash_mem(key, keylen);

  ient = hash % d->nalloc;

  de = d->list[ient];
  while(de) {
    if(de->hash == hash &&
       de->keylen == keylen &&
       !memcmp(de->key, key, keylen)) {
      if(value)
        *value = de->value;

      rv = 1;

      break;
    }

    de = de->next;
  }

  return(rv);
}

int dict_store (struct dict *d,
                void *key, uint32_t keylen,
                void *value) {
  uint32_t hash, ient, oent, new_nalloc, nadd;
  struct dict_entry **new_list = (struct dict_entry **) NULL;
  struct dict_entry *de = (struct dict_entry *) NULL;
  struct dict_entry *pde, *ode;

  hash = hash_mem(key, keylen);

  if(d->nused >= d->nalloc) {
    /* Extend, uses power of two algorithm to produce odd numbers
       (some of which are Mersenne primes). */
    if(d->nalloc)
      new_nalloc = d->nalloc * 2 + 1;
    else
      new_nalloc = 7;

    new_list
      = (struct dict_entry **) malloc(new_nalloc * sizeof(struct dict_entry *));
    if(!new_list)
      return(-1);

    for(ient = 0; ient < d->nalloc; ient++) {
      de = d->list[ient];
      while(de) {
        pde = de;
        de = de->next;
        pde->next = (struct dict_entry *) NULL;

        oent = pde->hash % new_nalloc;
        ode = new_list[oent];
        if(ode) {
          while(ode->next)
            ode = ode->next;

          ode->next = pde;
        }
        else
          new_list[oent] = pde;
      }
    }

    free((void *) d->list);

    d->list = new_list;
    d->nalloc = new_nalloc;
  }

  de = (struct dict_entry *) malloc(sizeof(struct dict_entry));
  if(!de)
    return(-1);

  nadd = 1;

  de->key = (void *) malloc(keylen);
  if(!de->key)
    return(-1);

  de->hash = hash;
  de->keylen = keylen;
  memcpy(de->key, key, keylen);
  de->value = value;
  de->next = (struct dict_entry *) NULL;

  ient = hash % d->nalloc;

  pde = d->list[ient];
  if(pde) {
    ode = (struct dict_entry *) NULL;
    for(;;) {
      if(pde->hash == hash &&
         pde->keylen == keylen &&
         !memcmp(pde->key, key, keylen)) {
        /* Replace */
        if(ode)
          ode->next = de;
        else
          d->list[ient] = de;

        de->next = pde->next;

        free((void *) pde->key);
        free((void *) pde);

        nadd = 0;

        break;
      }

      ode = pde;

      if(pde->next)
        pde = pde->next;
      else {
        pde->next = de;
        break;
      }
    }
  }
  else
    d->list[ient] = de;

  d->nused += nadd;

  return(0);
}

