#ifndef PATHSPACE_H
#define PATHSPACE_H

#include "corona_common.h"
#include "mf.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#ifndef PATHSPACE_MAX_VERTS
#define PATHSPACE_MAX_VERTS 32
// #define PATHSPACE_MAX_VERTS 3 // direct light only, remember to make clean.
#endif

// fix sampling dimensions for every vertex
typedef enum path_sample_dim_t
{
  // first vertex
  s_dim_lambda        = 2,
  s_dim_time          = 3,
  // from sensor
  s_dim_image_x       = 0,
  s_dim_image_y       = 1,
  s_dim_aperture_x    = 4,
  s_dim_aperture_y    = 5,
  s_dim_camid         = 6,  // choose camera for stereo renders
  s_dim_num_pt_beg    = 7,  // number of dimensions for starting at the sensor

  // from light
  s_dim_envmapvsarea  = 0,
  s_dim_lightsource   = 1,  // choose primitive
  s_dim_light_x       = 4,  // position on light
  s_dim_light_y       = 5,
  s_dim_edf_x         = 6,  // direction from light
  s_dim_edf_y         = 7,
  s_dim_num_lt_beg    = 8,  // number of dimensions for starting at the light

  // extend the path
  s_dim_free_path     = 0,  // travel in media (needs to be 0 because it's the only thing after sensor/emit)
  s_dim_scatter_mode  = 3,  // select lobe R/T
  s_dim_russian_r     = 4,  // russian roulette to terminate the path
  s_dim_omega_x       = 1,  // new direction
  s_dim_omega_y       = 2,  // 
  s_dim_num_extend    = 5,  // number of dims for path_extend()

  // next event estimation
  s_dim_nee_light1    = 0,
  s_dim_nee_light2    = 1,
  s_dim_nee_x         = 2,
  s_dim_nee_y         = 3,
  s_dim_num_nee       = 4,  // number of dims for next event
}
path_sample_dim_t;


// mode of scattering at a vertex
typedef enum vertex_scattermode_t
{
  s_reflect  = 1<<0,  // first two bits are mode of 
  s_transmit = 1<<1,  // reflection wrt normal
  s_volume   = 1<<2,  // volume scattering event (whole sphere, no cos)
  s_fiber    = 1<<3,  // fiber scattering (whole sphere, sin(fiber) instead of cos(normal))
  s_emit     = 1<<4,
  s_sensor   = 1<<5,

  s_diffuse  = 1<<6,  // then mode of bsdf wrt direction
  s_glossy   = 1<<7,
  s_specular = 1<<8,

  s_absorb   = 0,     // to signify none of the above
}
vertex_scattermode_t;


// store some stuff which is too small to justify a full int
typedef enum vertex_flags_t
{
  s_none        = 0,
  s_inside      = 1,  // only valid in sample() during construction in the same direction
  s_environment = 2,
}
vertex_flags_t;

typedef enum vertex_roughness_t
{
  s_rough_isotropic_beckmann = 0,
  s_rough_anisotropic_beckmann = 1,
}
vertex_roughness_t;

typedef struct vertex_shading_t
{
  // anisotropic roughness
  vertex_roughness_t roughness_type;
  float roughness, roughness_v, roughness_uv;
  mf_t rs;         // specular coefficient
  mf_t rd;         // diffuse coefficient
  mf_t rg;         // glossy coefficient
  mf_t em;         // emission
  float t;         // alpha-transparency
}
vertex_shading_t;

// encapsulates information for participating media
typedef struct vertex_volume_t
{
  mf_t mu_s;       // volume scattering coefficient
  mf_t mu_t;       // extinction coefficient 
  float mean_cos;  // mean cosine in medium
  mf_t ior;        // index of refraction
  int shader;      // volume shader or -1

  int num_lobes;   // number of extra lobes, if any
  mf_t l_mu_s[5];  // scattering
  mf_t l_mu_t[5];  // extinction
  float l_g[5];    // mean cosine/phase param
}
vertex_volume_t;

typedef enum vertex_manifold_type_t
{
  s_free = 0,
  s_pinned_position = 1,
  s_pinned_direction = 2,
}
vertex_manifold_type_t;

// stores stuff for differential geometry.
// has its name from wenzel's specular manifold exploration
typedef struct vertex_manifold_t
{
  // half vector to vertex area jacobian matrix, in sparse tridiagonal blocks
  // storage: the block matrices are 2x2 for surfaces and 3x3 for volumes. at
  // the boundaries, 2x3 and 3x2 are possible.
  // dh_i/dx_{i-1,i,i+1}, respectively.
  float a[9], b[9], c[9];
  float Tp[9]; // required for LU factorization to invert the full block-tridiagonal matrix
  float h[3]; // constraint: half vector in beckman [0] [1] or three components for volume scattering
  // then including volume constraint: distance( v[i-1], v[i] ) + distance ( v[i], v[i+1] )
  float dndu[3], dndv[3];
  float dpdu[3], dpdv[3];
  mf_t  eta;  // cached path_eta_ratio()
  vertex_manifold_type_t type;
}
vertex_manifold_t;

typedef struct vertex_t
{
  // geometric information
  hit_t hit;

  // monte carlo related
  mf_t pdf;               // on-surface probability (vertex area measure)
  mf_t pdf_mnee;          // XXX replace by generic pdf cache!
  mf_t throughput;        // throughput as constructed (incoming to this vertex)
  mf_t total_throughput;  // sum of throughputs up to when this vertex was the last one (includes emission at every previous, too)
  int rand_beg;           // first random number dimension used to construct this vertex
  int rand_cnt;           // number of random numbers drawn
  int tech;               // identifies technique used to create this vertex

  // slots for unified shading information
  vertex_shading_t shading;

  // volume properties of material inside the shape this vertex is on.
  // these will be set on the edge if transmitting inside.
  vertex_volume_t interior;

  // differential geometry data
  vertex_manifold_t diffgeo;

  vertex_flags_t flags;                 // various flags (envmap..)
  vertex_scattermode_t mode;            // reflect or transmit mode of /actually sampled or evaluated/ event.
  vertex_scattermode_t material_modes;  // potentially possible modes at this event, as far as material is concerned.
  vertex_scattermode_t culled_modes;    // light transport may choose to constrain material modes. affects pdf.
}
vertex_t;

// edge in a path, connects two vertices
typedef struct edge_t
{
  mf_t contribution;   // emissive volume contribution Le/p(v)
  mf_t transmittance;  // cached transmittance
  mf_t pdf;            // cached pdf to sample free path length
  float omega[3];      // direction pointing in the ray tracing direction from v[0] -> v[k]
  float dist;          // distance
  vertex_volume_t vol; // volume properties
}
edge_t;

// extra data to connect to a sensor.
// this belongs to the path because a path is only connected
// to one sensor. this way we keep the vertices all the same.
typedef struct sensor_t
{
  float pixel_i;      // pixel position
  float pixel_j;      //
  int   pixel_set;    // if != 0 force pixel_{i,j} which have been predetermined.
  float aperture_x;   // numbers used to create aperture sample.
  float aperture_y;   // these sometimes need to stay fixed for certain mlt mutations.
  int   aperture_set; // same as pixel_set but for aperture random numbers.
  int   camid;        // camera id to refer to list of cameras in view_t
  int   camid_set;    // keep already set camid
  float pixel_dir_i;  // store incoming direction on pixel for light field rendering
  float pixel_dir_j;
}
sensor_t;

// path struct bundling all the information of a transport path.
// be sure to change path_init() accordingly when messing with this!
typedef struct path_t
{
  mf_t lambda;       // wavelength
  mf_t throughput;   // path manipulation methods store the result for the final f/p here.
  int32_t length;    // number of vertices in this path
  float time;        // time, 0 is shutter open, 1 is shutter close
  float temperature; // used for tempering.
  uint64_t index;    // path index (for quasi-Monte Carlo numbers)
  float tangent_frame_scrambling; // random number which takes care that parametrisation flips in volume tangent frames are scrambled a bit
  sensor_t sensor;
  vertex_t v[PATHSPACE_MAX_VERTS];      // v[0] is start (light or sensor).
  edge_t   e[PATHSPACE_MAX_VERTS+1];    // e[i] leads up to v[i] (connecting from v[i-1].
  int debug_volume_bridge;
}
__attribute__((aligned(16))) // for faster memcpy (sse) that doesn't crash.
path_t;


// initialize everything to be an empty path.
void path_init(path_t *path, uint64_t index, int camid);

// copy only the relevant portions of memory (faster than assignment)
void path_copy(path_t *path, const path_t *other);

// sample new vertex. if length==0, sample new start point,
// from sensor or light, depending on v[0].mode & s_emit.
// the path is assumed to be initialized to all 0 (path_init)
// and will respect v[0].rand_beg (in case you do bidirectional)
// returns != 0 if the path sampling failed/the photon died
int path_extend(path_t *path);

// uses the russian roulette random number of the current last
// vertex of the path to determine whether the path should end
// here. that is, if this returns != 0, the path length should
// not be increased further. throughput and vertex pdf will be
// adjusted in both cases.
int path_russian_roulette(path_t *path, const float p_survival);

// pop last vertex of path, decreasing length by one.
// random number offsets will be maintained, so this action can be
// replayed in a kelemen metropolis context.
// only makes sense if there are more than 2 vertices.
void path_pop(path_t *path);

// return 1 if p->v[e] is visible from p->v[v-1] (along edge e)
int path_visible(path_t *p, int e);


// converts the outgoing throughput at v[v] to incoming throughput.
// also updates path contribution by taking account for point and segment emission.
// v[v].pdf is assumed to be in area measure.
void path_update_throughput(path_t* path, int v);

// get part of the measurement contribution without G (bsdfs times transmittances times W and L),
// that is, the measure space is projected solid angle.
// note that this is the expanded path integral (i.e. eq (8.7) p.223 in veach's thesis), not
// what veach defines as measurement contribution. the difference is that this accounts for
// light emission at every inner vertex of the path, too.
md_t path_measurement_contribution_dwp(path_t *path, int s, int e);

// get measurement contribution f in vertex area measure (this is how veach defines it)
// this calls path_measurement_contribution_dwp and multiplies the jacobian
// determinant |dwp/dx| (the geometric terms).
md_t path_measurement_contribution_dx(path_t *path, int s, int e);

// get pdf p as sampled, in product area measure
md_t path_pdf(const path_t *path);

// get throughput as sampled, X = f/p
mf_t path_throughput(const path_t *path);

// return on-surface pdf of vertex v if it had been sampled via path_extend at path->length = v
mf_t path_pdf_extend(const path_t *path, int v);

// return on-surface pdf of vertex v if it had been sampled the other way around via
// extension of the reverse path from v+1
mf_t path_pdf_extend_adjoint(const path_t *path, int v);

// reverse a path's vertices, make data consistent.
void path_reverse(path_t *path, const path_t *input);

// connect two paths, extending path1 by a connection edge and the reverse of path2.
mf_t path_connect(path_t *path1, const path_t *path2);

// deterministic path merging. biased. path1 will be altered, path1->v[path1->length-1] will stay.
// this is a deterministic, biased merge, if you're not doing this
// deterministically you should adjust the throughput by the pdf externally.
mf_t path_merge(path_t *path1, const path_t *path2);

typedef enum path_propagation_mode_t
{
  s_propagate_sample      = 0,  // regular free path sampling, as for vanilla path tracing (or light tracing)
  s_propagate_mutate      = 1,  // markov chain mutation mode, input pre-mutated distance (multichain mutation)
  s_propagate_reconstruct = 2,  // reconstruction mode, only verify exactly the path pointed to by the input and recreate all vertex data (photon maps)
}
path_propagation_mode_t;

// project changed vertex position at vertex v back to the surface/resample volume distances,
// by tracing a ray from vertex v-1
// returns 0 on success.
// if used in a metropolis context, pass mode=s_propagate_keep_dist to make it
// behave like a perturbation.
int path_project(path_t *path, int v, const path_propagation_mode_t mode);

// propagate a direction change in path->e[v] through the path
// by tracing a ray from vertex v-1 in direction e[v].omega
// returns 0 on success.
// if used in a metropolis context, pass mode=s_propagate_keep_dist to make it
// behave like a perturbation.
int path_propagate(path_t *path, int v, const path_propagation_mode_t mode);

// computes lambert's law for the given vertex and direction.
// this may seem trivial but has special cases for volumes and fibres.
float path_lambert(const path_t *p, int v, const float *omega);

// computes the standard geometric term of the given edge.
// cosines are included on demand (not for volumes) and the square falloff
// is not computed for environment map connections etc.
float path_G(const path_t *p, int e);

// initialise the volume properties on edge path->e[e].
// note that the vertex scattering mode path->v[e-1].mode needs to be set for
// this to work correctly (i.e. call shader_brdf() or sample before)
int path_edge_init_volume(path_t *path, int e);

// helper to find the ior ratio at the interface of vertex v.
// considers inside/outside, as well as nested volumes. will return n1/n2 where
// n1 is the ior of e[v] and n2 is the ior of e[v+1] if we transmit through the interface.
// only returns the correct thing after ior has been inited on v[v], which
// is /after/ the prepare shaders have been run (and the bsdf prepare in particular)
mf_t path_eta_ratio(const path_t *path, int v);

// pretty-print the path
void path_print(const path_t *path, FILE *f);

// pretty-print the worldspace to halfvector space matrix
void path_print_manifold_matrix(const path_t *path, FILE *f);

static inline void path_volume_vacuum(vertex_volume_t *v)
{
  memset(v, 0, sizeof(vertex_volume_t));
  v->shader = -1;
  v->ior = mf_set1(1.0f);
}

static inline void path_set_pixel(path_t *p, float i, float j)
{
  p->sensor.pixel_i = i;
  p->sensor.pixel_j = j;
  p->sensor.pixel_set = 1;
}

static inline void path_set_aperture(path_t *p, float u, float v)
{
  p->sensor.aperture_x = u;
  p->sensor.aperture_y = v;
  p->sensor.aperture_set = 1;
}

#endif
