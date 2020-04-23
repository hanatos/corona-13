/*
    This file is part of corona-13.
    copyright (c) 2015 johannes hanika

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "points.h"
#include <stdlib.h>

// xorshift128+, period 2^128-1, apparently passes all TestU01 suite tests.
typedef struct points_state_t
{
  uint64_t state0;
  uint64_t state1;
}
points_state_t;

typedef struct points_t
{
  points_state_t *s;
}
points_t;

points_t *points_init(const unsigned int num_threads, uint64_t frame)
{
  points_t *p = (points_t *)malloc(sizeof(points_t));
  p->s = (points_state_t *)malloc(sizeof(points_state_t)*num_threads);
  for(int k=0;k<num_threads;k++)
  {
    p->s[k].state0 = 1 + k;
    p->s[k].state1 = 2 + k + frame;
  }
  return p;
}

void points_cleanup(points_t *p)
{
  free(p->s);
  free(p);
}

void points_set_state(points_t *p, int thread_num, uint64_t s0, uint64_t s1)
{
  p->s[thread_num].state0 = 1 + thread_num + s0;
  p->s[thread_num].state1 = 2 + thread_num + s1;
  for(int k=0;k<10;k++) // discard a couple of rands to scramble state
    (void)points_rand(p, thread_num);
}

float points_rand(points_t *p, const unsigned int thread_num)
{
  uint64_t s1 = p->s[thread_num].state0;
  uint64_t s0 = p->s[thread_num].state1;
  p->s[thread_num].state0 = s0;
  s1 ^= s1 << 23;
  s1 ^= s1 >> 17;
  s1 ^= s0;
  s1 ^= s0 >> 26;
  p->s[thread_num].state1 = s1;
  // return (state0 + state1) / ((double)((uint64_t)-1) + 1.0);
  uint32_t v = 0x3f800000 | ((p->s[thread_num].state0+p->s[thread_num].state1)>>41); // faster than double version.
  return (*(float*)&v) - 1.0f;
}

void points_print_info(FILE *fd)
{
  fprintf(fd, "points   : xorshift128+\n");
}
