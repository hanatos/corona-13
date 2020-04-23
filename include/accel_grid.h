/*
    This file is part of corona-6: radiata.

    corona-6: radiata is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-6: radiata is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-6: radiata.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef GRID_H
#define GRID_H

#include <assert.h>
#ifdef __SSE__
  #include <xmmintrin.h>
#else
  #error "Grid needs SSE! turn it on in your compile flags."
#endif

// fisttpl ignores float rounding mode (=> (int))
#define FLOAT_TO_INT(in,out)  \
                    __asm__ __volatile__ ("fistpl %0" : "=m" (out) : "t" (in) : "st") ;

typedef struct hashcell_t
{
  int prim;
  int num_prims;
}
hashcell_t;

typedef struct accel_t
{
  int num_prims;
  int num_entries;
  int num_cells[3];
  float grid_size;
  float grid_factor;
  int index_size;
  int *index;
  prim_t *prim;
  hashcell_t *hashtable;
  unsigned int *gridmask;
  float aabb[6];

  int cachehits;
  int cachelookups;
}
accel_t;

#define grid_t accel_t

static inline unsigned int accel_populated(const grid_t *g, const int *i)
{
  const unsigned int index = g->num_cells[0]*(g->num_cells[1]*i[2] + i[1]) + i[0];
  return g->gridmask[index>>5] & (1<<(index&31));
}

static inline void accel_populate(grid_t *g, const int *i)
{
  assert(i[0] < g->num_cells[0] && i[0] >= 0);
  assert(i[1] < g->num_cells[1] && i[1] >= 0);
  assert(i[2] < g->num_cells[2] && i[2] >= 0);
  const unsigned int index = g->num_cells[0]*(g->num_cells[1]*i[2] + i[1]) + i[0];
  g->gridmask[index>>5] |= (1<<(index&31));
}

static inline void accel_gridpos(const grid_t *g, const float *x, int *i)
{
  //for(int k=0;k<3;k++) i[k] = (int)((x[k] - g->aabb[k])*g->grid_factor);
  for(int k=0;k<3;k++) FLOAT_TO_INT((x[k] - g->aabb[k])*g->grid_factor, i[k]);
}

#if 0
typedef  unsigned long  int  ub4;   /* unsigned 4-byte quantities */
typedef  unsigned       char ub1;   /* unsigned 1-byte quantities */

#define hashsize(n) ((ub4)1<<(n))
#define hashmask(n) (hashsize(n)-1)

/*
--------------------------------------------------------------------
mix -- mix 3 32-bit values reversibly.
For every delta with one or two bits set, and the deltas of all three
  high bits or all three low bits, whether the original value of a,b,c
  is almost all zero or is uniformly distributed,
* If mix() is run forward or backward, at least 32 bits in a,b,c
  have at least 1/4 probability of changing.
* If mix() is run forward, every bit of c will change between 1/3 and
  2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
mix() was built out of 36 single-cycle latency instructions in a 
  structure that could supported 2x parallelism, like so:
      a -= b; 
      a -= c; x = (c>>13);
      b -= c; a ^= x;
      b -= a; x = (a<<8);
      c -= a; b ^= x;
      c -= b; x = (b>>13);
      ...
  Unfortunately, superscalar Pentiums and Sparcs can't take advantage 
  of that parallelism.  They've also turned some of those single-cycle
  latency instructions into multi-cycle latency instructions.  Still,
  this is the fastest good hash I could find.  There were about 2^^68
  to choose from.  I only looked at a billion or so.
--------------------------------------------------------------------
*/
#define mix(a,b,c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

/*
--------------------------------------------------------------------
hash() -- hash a variable-length key into a 32-bit value
  k       : the key (the unaligned variable-length array of bytes)
  len     : the length of the key, counting by bytes
  initval : can be any 4-byte value
Returns a 32-bit value.  Every bit of the key affects every bit of
the return value.  Every 1-bit and 2-bit delta achieves avalanche.
About 6*len+35 instructions.

The best hash table sizes are powers of 2.  There is no need to do
mod a prime (mod is sooo slow!).  If you need less than 32 bits,
use a bitmask.  For example, if you need only 10 bits, do
  h = (h & hashmask(10));
In which case, the hash table should have hashsize(10) elements.

If you are hashing n strings (ub1 **)k, do it like this:
  for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);

By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
code any way you wish, private, educational, or commercial.  It's free.

See http://burtleburtle.net/bob/hash/evahash.html
Use for hash table lookup, or anything where one collision in 2^^32 is
acceptable.  Do NOT use for cryptographic purposes.
--------------------------------------------------------------------
*/

ub4 hash( k, length, initval)
register ub1 *k;        /* the key */
register ub4  length;   /* the length of the key */
register ub4  initval;  /* the previous hash, or an arbitrary value */
{
   register ub4 a,b,c,len;

   /* Set up the internal state */
   len = length;
   a = b = 0x9e3779b9;  /* the golden ratio; an arbitrary value */
   c = initval;         /* the previous hash value */

   /*---------------------------------------- handle most of the key */
   while (len >= 12)
   {
      a += (k[0] +((ub4)k[1]<<8) +((ub4)k[2]<<16) +((ub4)k[3]<<24));
      b += (k[4] +((ub4)k[5]<<8) +((ub4)k[6]<<16) +((ub4)k[7]<<24));
      c += (k[8] +((ub4)k[9]<<8) +((ub4)k[10]<<16)+((ub4)k[11]<<24));
      mix(a,b,c);
      k += 12; len -= 12;
   }

   /*------------------------------------- handle the last 11 bytes */
   c += length;
   switch(len)              /* all the case statements fall through */
   {
   case 11: c+=((ub4)k[10]<<24);
   case 10: c+=((ub4)k[9]<<16);
   case 9 : c+=((ub4)k[8]<<8);
      /* the first byte of c is reserved for the length */
   case 8 : b+=((ub4)k[7]<<24);
   case 7 : b+=((ub4)k[6]<<16);
   case 6 : b+=((ub4)k[5]<<8);
   case 5 : b+=k[4];
   case 4 : a+=((ub4)k[3]<<24);
   case 3 : a+=((ub4)k[2]<<16);
   case 2 : a+=((ub4)k[1]<<8);
   case 1 : a+=k[0];
     /* case 0: nothing left to add */
   }
   mix(a,b,c);
   /*-------------------------------------------- report the result */
   return c;
}
#endif

static inline hashcell_t *accel_hashcelli(const grid_t *g, const int *i)
{
  //return g->hashtable + ((3*i[0] + 5*i[1] + 7*i[2]) & (g->num_entries - 1));
  //return g->hashtable + (((i[0]*73856093) ^ (i[1]*19349663) ^ (i[2]*83492791)) & (g->num_entries - 1));//% g->num_entries;
  //return g->hashtable + ((((i[0]<<17) ^ (i[1]<<7) ^ (i[2]<<3))) & (g->num_entries - 1));//% g->num_entries;
  return g->hashtable + ((g->num_cells[0]*(g->num_cells[1]*i[2] + i[1]) + i[0]) & (g->num_entries - 1));
  // sort-of-MMH:
  /*unsigned long long int sum = 0;
  const int key = 12345;
  for(int k=0;k<3;k++)
  {
    sum += i[k] * key;
  }
  sum = sum % (unsigned long long int)((g->num_entries - 1) + 15ULL);
  sum &= (g->num_entries - 1);
  return g->hashtable + sum;*/
  // bob
  //return g->hashtable + (hash((ub1*)i, 9*sizeof(float), 1234) & (g->num_entries - 1));
}

static inline hashcell_t *accel_hashcell(grid_t *g, const float *x)
{
  int i[3];
  accel_gridpos(g, x, i);
  return accel_hashcelli(g, i);
}

accel_t* accel_init(const int num_prims, prim_t *prim);
void accel_cleanup(struct accel_t *g);
void accel_build(struct accel_t *g, const char* filename);
void accel_intersect(const struct accel_t *g, const ray_t *ray, rayhit_t *hit);
int  accel_visible(const struct accel_t *g, const ray_t *ray);

static inline void accel_print_info(FILE *fd)
{
  fprintf(fd, "accel    : regular grid + hashtable\n");
  fprintf(fd, "           hashtable has %d entries\n", rt.accel->num_entries);
  fprintf(fd, "           %dx%dx%d grid\n", rt.accel->num_cells[0], rt.accel->num_cells[1], rt.accel->num_cells[2]);
  fprintf(fd, "           grid size: %f \n", rt.accel->grid_size);
}


#endif
