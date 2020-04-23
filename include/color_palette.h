#pragma once

static float rye[5][3] = {
  { 0.0f, 0.0f, 0.5f },
  { 0.0f, 1.0f, 0.0f },
  { 1.0f, 1.0f, 0.0f },
  { 1.0f, 0.5f, 0.0f },
  { 1.0f, 0.0f, 0.0f },
};

static inline void lerp(
    float w, 
    const float c0[3], 
    const float c1[3], 
    float c[3])
{
  c[0] = w * c1[0] + (1.f - w) * c0[0];
  c[1] = w * c1[1] + (1.f - w) * c0[1];
  c[2] = w * c1[2] + (1.f - w) * c0[2];
}

static inline void color_palette(float w, float color[3])
{
  if (w < 0.25f) 
  {
    w /= 0.25;
    lerp(w, rye[0], rye[1], color);
  }
  else if (w < 0.5f)
  {
    w -= 0.25;
    w /= 0.25;
    lerp(w, rye[1], rye[2], color);
  }
  else if (w < 0.75)
  {
    w -= 0.5;
    w /= 0.25;
    lerp(w, rye[2], rye[3], color);
  }
  else 
  {
    w -= 0.75;
    w /= 0.25;
    lerp(w, rye[3], rye[4], color);
  }
}
