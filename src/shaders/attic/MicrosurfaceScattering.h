#ifndef _MICROSURFACESCATTERING_
#define _MICROSURFACESCATTERING_

#include <glm/glm.hpp>
using namespace glm;
#include <complex>


/************* MICROSURFACE HEIGHT DISTRIBUTION *************/

/* API */
class MicrosurfaceHeight
{
public:
	// height PDF	
	virtual float P1(const float h) const=0; 
	// height CDF	
	virtual float C1(const float h) const=0; 
	// inverse of the height CDF
	virtual float invC1(const float U) const=0; 
};

/* Uniform height distribution in [-1, 1] */
class MicrosurfaceHeightUniform : public MicrosurfaceHeight
{
public:	
	// height PDF	
	virtual float P1(const float h) const; 
	// height CDF	
	virtual float C1(const float h) const; 
	// inverse of the height CDF
	virtual float invC1(const float U) const; 
};

/* Gaussian height distribution N(0,1) */
class MicrosurfaceHeightGaussian : public MicrosurfaceHeight
{
public:	
	// height PDF	
	virtual float P1(const float h) const; 
	// height CDF	
	virtual float C1(const float h) const; 
	// inverse of the height CDF
	virtual float invC1(const float U) const; 
};


/************* MICROSURFACE SLOPE DISTRIBUTION *************/

/* API */
class MicrosurfaceSlope
{
public:
	MicrosurfaceSlope(const float alpha_x=1.0f, const float alpha_y=1.0f)
		: m_alpha_x(alpha_x), m_alpha_y(alpha_y)
	{}

public:
	// roughness
	const float m_alpha_x, m_alpha_y;
	// projected roughness in wi
	float alpha_i(const vec3& wi) const; 

public:
	// distribution of normals (NDF)	
	float D(const vec3& wm) const; 
	// distribution of visible normals (VNDF)
	float D_wi(const vec3& wi, const vec3& wm) const; 
	// sample the VNDF
	vec3 sampleD_wi(const vec3& wi, const float U1, const float U2) const;

public:
	// distribution of slopes
	virtual float P22(const float slope_x, const float slope_y) const=0; 
	// Smith's Lambda function
	virtual float Lambda(const vec3& wi) const=0;
	// projected area towards incident direction
	virtual float projectedArea(const vec3& wi) const=0;
	// sample the distribution of visible slopes with alpha=1.0
	virtual vec2 sampleP22_11(const float theta_i, const float U1, const float U2) const=0;
};

/* Beckmann slope distribution */
class MicrosurfaceSlopeBeckmann : public MicrosurfaceSlope
{
public:
	MicrosurfaceSlopeBeckmann(const float alpha_x=1.0f, const float alpha_y=1.0f)
		: MicrosurfaceSlope(alpha_x, alpha_y)
	{}

	// distribution of slopes
	virtual float P22(const float slope_x, const float slope_y) const; 
	// Smith's Lambda function
	virtual float Lambda(const vec3& wi) const;
	// projected area towards incident direction
	virtual float projectedArea(const vec3& wi) const;
	// sample the distribution of visible slopes with alpha=1.0
	virtual vec2 sampleP22_11(const float theta_i, const float U1, const float U2) const;
};

/* GGX slope distribution */
class MicrosurfaceSlopeGGX : public MicrosurfaceSlope
{
public:
	MicrosurfaceSlopeGGX(const float alpha_x=1.0f, const float alpha_y=1.0f)
		: MicrosurfaceSlope(alpha_x, alpha_y)
	{}

	// distribution of slopes
	virtual float P22(const float slope_x, const float slope_y) const; 
	// Smith's Lambda function
	virtual float Lambda(const vec3& wi) const;
	// projected area towards incident direction
	virtual float projectedArea(const vec3& wi) const;
	// sample the distribution of visible slopes with alpha=1.0
	virtual vec2 sampleP22_11(const float theta_i, const float U1, const float U2) const;
};


/************* MICROSURFACE *************/

/* API */
class Microsurface 
{
public:
	// height distribution
	const MicrosurfaceHeight* m_microsurfaceheight; 
	// slope distribution
	const MicrosurfaceSlope* m_microsurfaceslope; 
  const MicrosurfaceHeightUniform m_h;
  const MicrosurfaceSlopeGGX m_s;

public:

	Microsurface(const bool height_uniform, // uniform or Gaussian height distribution
				const bool slope_beckmann, // Beckmann or GGX slope distribution
				const float alpha_x,
				const float alpha_y) : 
    m_h(),
    m_s(alpha_x, alpha_y)
		// m_microsurfaceheight((height_uniform) ? 
		//   static_cast<MicrosurfaceHeight*>(new MicrosurfaceHeightUniform) 
		// : static_cast<MicrosurfaceHeight*>(new MicrosurfaceHeightGaussian)),
		// m_microsurfaceslope((slope_beckmann) ? 
		//   static_cast<MicrosurfaceSlope*>(new MicrosurfaceSlopeBeckmann(alpha_x, alpha_y)) 
		// : static_cast<MicrosurfaceSlope*>(new MicrosurfaceSlopeGGX(alpha_x, alpha_y)))
	{
    m_microsurfaceheight = &m_h;
    m_microsurfaceslope = &m_s;
  }
		
	~Microsurface() 
	{
		// delete m_microsurfaceheight; 
		// delete m_microsurfaceslope;
	}

	// evaluate BSDF with a random walk (stochastic but unbiased)
	// scatteringOrder=0 --> contribution from all scattering events
	// scatteringOrder=1 --> contribution from 1st bounce only
	// scatteringOrder=2 --> contribution from 2nd bounce only, etc..
	virtual float eval(const vec3& wi, const vec3& wo, const int scatteringOrder=0) const; 

	// sample BSDF with a random walk
	// scatteringOrder is set to the number of bounces computed for this sample
	virtual float sample(const vec3& wi, vec3 &wo, int& scatteringOrder, float r0, float r1, float r2) const;
	// vec3 sample(const vec3& wi) const {int scatteringOrder; return sample(wi, scatteringOrder);}
  virtual float pdf(const vec3& wi, const vec3& wo) const;

public:
	// masking function
	float G_1(const vec3& wi) const;
	// masking function at height h0
	float G_1(const vec3& wi, const float h0) const;
	// sample height in outgoing direction
	float sampleHeight(const vec3& wo, const float h0, const float U) const;
	float sampleHeight2(const vec3& wo, const float h0, const float U, float &thr) const;

public:
	// evaluate local phase function 
	virtual float evalPhaseFunction(const vec3& wi, const vec3& wo) const=0;
	// sample local phase function
  virtual float samplePhaseFunction(const vec3& wi, vec3 &wo, float U0, float U1, float U2) const = 0;

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	virtual float evalSingleScattering(const vec3& wi, const vec3& wo) const=0; 
};

/* Microsurface made of conductor material */
class MicrosurfaceConductor : public Microsurface
{
  std::complex<float> m_ior;
public:
	MicrosurfaceConductor(const bool height_uniform, // uniform or Gaussian
				 const bool slope_beckmann, // Beckmann or GGX
				 const float alpha_x,
				 const float alpha_y,
         const std::complex<float> ior)
		: Microsurface(height_uniform, slope_beckmann, alpha_x, alpha_y)
    , m_ior(ior)
	{}

public:
	// evaluate local phase function 
	float evalPhaseFunction(const vec3& wi, const vec3& wo) const;
	// sample local phase function
  float samplePhaseFunction(const vec3& wi, vec3 &wo, float U0, float U1, float U2) const;

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	float evalSingleScattering(const vec3& wi, const vec3& wo) const; 
};

/* Microsurface made of conductor material */
class MicrosurfaceDielectric : public Microsurface
{
public:
	const float m_eta;
public:
	MicrosurfaceDielectric(const bool height_uniform, // uniform or Gaussian
				 const bool slope_beckmann, // Beckmann or GGX
				 const float alpha_x,
				 const float alpha_y,
				 const float eta = 1.5f)
		: Microsurface(height_uniform, slope_beckmann, alpha_x, alpha_y),
		m_eta(eta)
	{}

	// evaluate BSDF with a random walk (stochastic but unbiased)
	// scatteringOrder=0 --> contribution from all scattering events
	// scatteringOrder=1 --> contribution from 1st bounce only
	// scatteringOrder=2 --> contribution from 2nd bounce only, etc..
	float eval(const vec3& wi, const vec3& wo, const int scatteringOrder=0) const; 

	// sample final BSDF with a random walk
	// scatteringOrder is set to the number of bounces computed for this sample
	float sample(const vec3& wi, vec3& wo, int& scatteringOrder, float r0, float r1, float r2) const;
  float pdf(const vec3& wi, const vec3& wo) const;

public:
	// evaluate local phase function 
	float evalPhaseFunction(const vec3& wi, const vec3& wo) const;
	float evalPhaseFunction(const vec3& wi, const vec3& wo, const bool wi_outside, const bool wo_outside) const;
  // not used:
  float samplePhaseFunction(const vec3& wi, vec3 &wo, float U0, float U1, float U2) const {return 0;}
	// sample local phase function
	float samplePhaseFunction(const vec3& wi, vec3 &wo, const bool wi_outside, bool& wo_outside, float U0, float U1, float U2) const; 

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	float evalSingleScattering(const vec3& wi, const vec3& wo) const; 

protected:
	float Fresnel(const vec3& wi, const vec3& wm, const float eta) const;
};

/* Microsurface made of conductor material */
class MicrosurfaceDiffuse : public Microsurface
{
  float m_reflectivity;
public:
	MicrosurfaceDiffuse(const bool height_uniform, // uniform or Gaussian
				 const bool slope_beckmann, // Beckmann or GGX
				 const float alpha_x,
				 const float alpha_y,
         const float reflectivity)
		: Microsurface(height_uniform, slope_beckmann, alpha_x, alpha_y)
    , m_reflectivity(reflectivity)
	{}

public:
	// evaluate local phase function 
	float evalPhaseFunction(const vec3& wi, const vec3& wo) const;
	// sample local phase function
	float samplePhaseFunction(const vec3& wi, vec3& wo, float U0, float U1, float U2) const;

	// evaluate BSDF limited to single scattering 
	// this is in average equivalent to eval(wi, wo, 1);
	float evalSingleScattering(const vec3& wi, const vec3& wo) const; 
};


#endif

