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

#include "points.h"

typedef struct points_t { int none; } points_t;
static unsigned int points_seedvalue = 0x15A4E35;
points_t *points_init(const unsigned int num_threads, uint64_t frame) { points_seedvalue = (unsigned int)frame; return 0; }
void points_cleanup(points_t *p) { }

float points_rand(points_t *p, const unsigned int thread_num)
{
#pragma omp atomic
  points_seedvalue *= (unsigned int) 0x15A4E35;
  return (float) (points_seedvalue >> 16) / (0xFFFF + 1.0f);
}

void points_print_info(FILE *fd)
{
  fprintf(fd, "points   : 1337 rand\n");
}

void points_set_state(points_t *p, int thread_num, uint64_t s0, uint64_t s1) {}
