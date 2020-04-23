// this is somewhat ugly, but abstracts the measurement contribution computation in different measurement spaces nicely.
// if you include this file in a function body as it is, it will compute the accumulated measurement contribution for
// projected solid angle measure. these two macros are entry points to compute per-vertex and per-edge jacobians to
// transform inte different spaces:
#ifndef PER_EDGE_JACOBIAN
#define PER_EDGE_JACOBIAN
#endif
#ifndef PER_VERTEX_JACOBIAN
#define PER_VERTEX_JACOBIAN
#endif
// first and last vertices of sub-path to eval measurement on:
#if !defined(MEASUREMENT_BEG) && !defined(MEASUREMENT_END)
  // if those aren't set, we check whether the path is complete at all:
  if(path->length < 2) return 0.0;
  if(!((path->v[0].mode & s_sensor) || (path->v[path->length-1].mode & s_sensor)))
    return 0.0;
#endif
#ifndef MEASUREMENT_BEG
#define MEASUREMENT_BEG 0
#endif
#ifndef MEASUREMENT_END
#define MEASUREMENT_END (path->length-1)
#endif

  md_t f = md_set1(1.0);
  md_t result = md_set1(0.0f);
  // loop goes over all edge indices from sensor to light
  const int lt = path->v[0].mode & s_emit;
  const int beg = lt ? MEASUREMENT_END : MEASUREMENT_BEG+1;
  const int end = lt ? MEASUREMENT_BEG : MEASUREMENT_END+1;
  const int inc = lt ? -1 : 1;
  for(int e=beg;e!=end;e+=inc)
  {
    const int v = lt ? e-1 : e; // vertex after this edge (closer to light)

    PER_EDGE_JACOBIAN // macro to mul stuff to f

    // volume attenuation
    if(mf_any(mf_gt(path->e[e].transmittance, mf_set1(0.0f))))
      f = md_mul(f, mf_2d(path->e[e].transmittance));
    else
      f = md_mul(f, mf_2d(shader_vol_transmittance(path, e)));

    // only inner vertices:
    if(e != end-inc)
    {
      PER_VERTEX_JACOBIAN // macro to mul stuff to f
      result = md_add(result, md_mul(f, mf_2d(lights_eval_vertex(path, v))));
      f = md_mul(f, mf_2d(shader_brdf(path, v)));
    }
    if(mf_all(mf_lte(md_2f(f), mf_set1(0.0f)))) break;
  }

  // eval light source L at last vertex
  if(path->v[MEASUREMENT_BEG].mode & s_emit)
    f = md_mul(f, mf_2d(lights_eval_vertex(path, MEASUREMENT_BEG)));
  if(path->v[MEASUREMENT_END].mode & s_emit)
    f = md_mul(f, mf_2d(lights_eval_vertex(path, MEASUREMENT_END)));

  // don't add to result in case path was absorbed
  if(path->v[MEASUREMENT_BEG].mode != s_absorb &&
     path->v[MEASUREMENT_END].mode != s_absorb)
    result = md_add(result, f);

  // eval camera W, shared for all path lengths:
  if((path->v[MEASUREMENT_BEG].mode & s_sensor) ||
     (path->v[MEASUREMENT_END].mode & s_sensor))
    result = md_mul(result, mf_2d(view_cam_eval(path)));
  return result;
#undef PER_EDGE_JACOBIAN
#undef PER_VERTEX_JACOBIAN
#undef MEASUREMENT_BEG
#undef MEASUREMENT_END
