/*
    This file is part of darktable,
    copyright (c) 2012 johannes hanika.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <assert.h>
#include <float.h>

// simple implementation of a heap/priority queue, using uint64_t as key and
// float values to sort the elements.
// meant to support scheduling of background jobs with priorities.
typedef struct heap_t
{
  uint64_t size;
  uint64_t end;
  uint64_t *keys;
  float    *vals;
  int      fixed_numbers; // if this is set, keys will be numbered from 1..size and immutable (useful for custom memory pools on the side)
}
heap_t;

static inline heap_t *heap_init(uint64_t size)
{
  heap_t *h = (heap_t *)malloc(sizeof(heap_t));
  h->keys = (uint64_t *)malloc(sizeof(uint64_t)*size);
  h->vals = (float    *)malloc(sizeof(float)   *size);
  h->size = size;
  h->end = 0;
  h->fixed_numbers = 0;
  return h;
}

static inline void heap_init_fixed(heap_t *h)
{
  h->fixed_numbers = 1;
  for(uint64_t k=0;k<h->size;k++) h->keys[k] = k;
}

static inline void heap_cleanup(heap_t *h)
{
  free(h->keys);
  free(h->vals);
  free(h);
}

static inline int heap_empty(heap_t *h)
{
  return h->end == 0;
}

static inline void heap_clear(heap_t *h)
{
  h->end = 0;
}

static inline int heap_full(heap_t *h)
{
  return h->end >= h->size;
}

static inline uint64_t _heap_parent(uint64_t i)
{
  return (i-1)/2;
}

static inline uint64_t _heap_child(uint64_t i, uint64_t right)
{
  return 2*i + 1 + right;
}

static inline void _heap_swap(heap_t *h, uint64_t i, uint64_t j)
{
  uint64_t tmpi = h->keys[i];
  h->keys[i] = h->keys[j];
  h->keys[j] = tmpi;

  float tmpf = h->vals[i];
  h->vals[i] = h->vals[j];
  h->vals[j] = tmpf;
}

// insert new element, growing the heap by one
// key is unused in case of fixed initialisation, and the correct one will be returned
static inline uint64_t heap_insert(heap_t *h, uint64_t key, float val)
{
  uint64_t pos = (h->end)++;
  assert(pos < h->size);
  if(!h->fixed_numbers)
    h->keys[pos] = key;
  else key = h->keys[pos];
  h->vals[pos] = val;

  while(pos >= 1)
  {
    uint64_t prt = _heap_parent(pos);
    if(h->vals[prt] < h->vals[pos])
    {
      _heap_swap(h, prt, pos);
      pos = prt;
    }
    else
    {
      break;
    }
  }
  return key;
}

// remove largest element from the heap
static inline void heap_remove(heap_t *h, uint64_t *key, float *val)
{
  *key = h->keys[0];
  *val = h->vals[0];

  if(--(h->end) == 0) return;
  _heap_swap(h, 0, h->end);

  uint64_t pos = 0;
  while(1)
  {
    uint64_t largest = pos;
    uint64_t c0 = _heap_child(pos, 0);
    uint64_t c1 = _heap_child(pos, 1);
    if(c0 < h->end && h->vals[c0] > h->vals[pos]) largest = c0;
    if(c1 < h->end && h->vals[c1] > h->vals[largest]) largest = c1;
    if(largest != pos)
    {
      _heap_swap(h, largest, pos);
      pos = largest;
    }
    else break;
  }
}

// remove largest element from the heap
static inline float heap_max_val(heap_t *h)
{
  if (h->end == 0)
    return -FLT_MAX;
  return h->vals[0];
}
