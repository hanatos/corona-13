/*
    This file is part of corona-13.

    corona-13 is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    corona-13 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with corona-13.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "pathspace.h"

// forward declare struct
struct sampler_t;
typedef struct sampler_t sampler_t;

sampler_t *sampler_init();
void sampler_cleanup(sampler_t *s);
void sampler_prepare_frame(sampler_t *s);
void sampler_clear(sampler_t *s);

// create a group of paths. when a valid connection between sensor and lights
// has been found, this shall call pointsampler_splat()
void sampler_create_path(path_t *path);

// this callback is needed to work with veach mlt large steps.
// it pretends to have sampled the given path with this method,
// and computes the throughput. for numerical reasons we use
// throughput = measurement contrib/product vertex area pdf directly,
// and not the pdf function below.
mf_t sampler_throughput(path_t *path);

// returns the multiple importance sampling weight of a path.
// this makes use of the technique (as sampled) stored at every vertex
// to find out with which technique the path should have been constructed,
// as well as knowledge about what other techniques this sampler
// combines, and the mis heuristic (mostly power i suppose). mostly this will be:
// pdf_t(path)^2 / sum_i { pdf_i(path)^2 }
md_t sampler_mis_weight(path_t *path);

// return the sum of all possible construction pdf, in projected solid angle measure.
// that is, for pt+nee, it would be the sum of extend + nee. to be used in the balance
// heuristic MIS context.
md_t sampler_sum_pdf_dwp(path_t *p);

// output info for sidecar file
void sampler_print_info(FILE *fd);
