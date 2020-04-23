#include "MicrosurfaceScattering.h"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <random>
extern "C"{
#include "corona_common.h"
#include "points.h"
}
using namespace std;

#ifndef M_PI
#define M_PI			3.14159265358979323846f	/* pi */
#define M_PI_2			1.57079632679489661923f	/* pi/2 */
#endif
#define INV_M_PI		0.31830988618379067153f /* 1/pi */
#define SQRT_M_PI		1.77245385090551602729f /* sqrt(pi) */
#define SQRT_2			1.41421356237309504880f /* sqrt(2) */
#define INV_SQRT_M_PI	0.56418958354775628694f /* 1/sqrt(pi) */
#define INV_2_SQRT_M_PI	0.28209479177387814347f /* 0.5/sqrt(pi) */
#define INV_SQRT_2_M_PI 0.3989422804014326779f /* 1/sqrt(2*pi) */
#define INV_SQRT_2		0.7071067811865475244f /* 1/sqrt(2) */


static bool IsFiniteNumber(float x) 
{
	return (x <= FLT_MAX && x >= -FLT_MAX); 
} 


// static std::random_device rd;
// static std::mt19937 gen(rd());
// static std::uniform_real_distribution<> dis(0, 1);

static float generateRandomNumber()
{
  return points_get_comp(rt.points, common_get_threadid());
	// return (float)dis(gen);
}


#if 0
static double erf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
 
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);
 
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    return sign*y;
}
#endif

static float erfinv(float x)
{
float w, p;
w = - logf((1.0f-x)*(1.0f+x));
if ( w < 5.000000f ) {
w = w - 2.500000f;
p = 2.81022636e-08f;
p = 3.43273939e-07f + p*w;
p = -3.5233877e-06f + p*w;
p = -4.39150654e-06f + p*w;
p = 0.00021858087f + p*w;
p = -0.00125372503f + p*w;
p = -0.00417768164f + p*w;
p = 0.246640727f + p*w;
p = 1.50140941f + p*w;
}
else {
w = sqrtf(w) - 3.000000f;
p = -0.000200214257f;
p = 0.000100950558f + p*w;
p = 0.00134934322f + p*w;
p = -0.00367342844f + p*w;
p = 0.00573950773f + p*w;
p = -0.0076224613f + p*w;
p = 0.00943887047f + p*w;
p = 1.00167406f + p*w;
p = 2.83297682f + p*w;
}
return p*x;
}



#ifndef M_PI
#define M_PI			3.14159265358979323846f	/* pi */
#endif

/* 
 * A method to compute the gamma() function.
 *
 */  

static double  abgam (double x)
{
  double  gam[10],
          temp;

  gam[0] = 1./ 12.;
  gam[1] = 1./ 30.;
  gam[2] = 53./ 210.;
  gam[3] = 195./ 371.;
  gam[4] = 22999./ 22737.;
  gam[5] = 29944523./ 19733142.;
  gam[6] = 109535241009./ 48264275462.;
  temp = 0.5*log (2*M_PI) - x + (x - 0.5)*log (x)
    + gam[0]/(x + gam[1]/(x + gam[2]/(x + gam[3]/(x + gam[4] /
	  (x + gam[5]/(x + gam[6]/x))))));

  return temp;
}

static double  ggamma (double x)
{
  double  result;
  result = exp (abgam (x + 5))/(x*(x + 1)*(x + 2)*(x + 3)*(x + 4));
  return result;
}

static double  beta (double m, double n)
{
  return (ggamma (m)*ggamma (n)/ggamma (m + n));
}


/************* MICROSURFACE HEIGHT DISTRIBUTION *************/

float MicrosurfaceHeightUniform::P1(const float h) const
{
	const float value = (h >= -1.0f && h <= 1.0f) ? 0.5f : 0.0f;
	return value;
}

float MicrosurfaceHeightUniform::C1(const float h) const
{
	const float value = std::min(1.0f, std::max(0.0f, 0.5f*(h+1.0f)));
	return value;
}

float MicrosurfaceHeightUniform::invC1(const float U) const
{
	const float h = std::max(-1.0f, std::min(1.0f, 2.0f*U-1.0f));
	return h;	
}

float MicrosurfaceHeightGaussian::P1(const float h) const
{
	const float value = INV_SQRT_2_M_PI * expf(-0.5f * h*h);
	return value;
}

float MicrosurfaceHeightGaussian::C1(const float h) const
{
	const float value = 0.5f + 0.5f * (float)erf(INV_SQRT_2*h);
	return value;
}

float MicrosurfaceHeightGaussian::invC1(const float U) const
{
	const float h = SQRT_2 * erfinv(2.0f*U - 1.0f);
	return h;	
}


/************* MICROSURFACE SLOPE DISTRIBUTION *************/

float MicrosurfaceSlope::D(const vec3& wm) const {
	if( wm.z <= 0.0f)
		return 0.0f;

	// slope of wm
	const float slope_x = -wm.x/wm.z;
	const float slope_y = -wm.y/wm.z;

	// value
	const float value = P22(slope_x, slope_y) / (wm.z*wm.z*wm.z*wm.z);
	return value;
}

float MicrosurfaceSlope::D_wi(const vec3& wi, const vec3& wm) const {
	if( wm.z <= 0.0f)
		return 0.0f;

	// normalization coefficient
	const float projectedarea = projectedArea(wi);
	if(projectedarea == 0)
		return 0;
	const float c = 1.0f / projectedarea;

	// value
	const float value = c * std::max(0.0f, dot(wi, wm)) * D(wm);
	return value;
}

vec3 MicrosurfaceSlope::sampleD_wi(const vec3& wi, const float U1, const float U2) const {

	// stretch to match configuration with alpha=1.0	
	const vec3 wi_11 = normalize(vec3(m_alpha_x * wi.x, m_alpha_y * wi.y, wi.z));

	// sample visible slope with alpha=1.0
	vec2 slope_11 = sampleP22_11(acosf(wi_11.z), U1, U2);

	// align with view direction
	const float phi = atan2(wi_11.y, wi_11.x);
	vec2 slope(cosf(phi)*slope_11.x - sinf(phi)*slope_11.y, sinf(phi)*slope_11.x + cos(phi)*slope_11.y);

	// stretch back
	slope.x *= m_alpha_x;
	slope.y *= m_alpha_y;

	// compute normal
	const vec3 wm = normalize(vec3(-slope.x, -slope.y, 1.0f));
	return wm;
}

float MicrosurfaceSlope::alpha_i(const vec3& wi) const
{
	const float invSinTheta2 = 1.0f / (1.0f - wi.z*wi.z);
	const float cosPhi2 = wi.x*wi.x*invSinTheta2;
	const float sinPhi2 = wi.y*wi.y*invSinTheta2;
	const float alpha_i = sqrtf( cosPhi2*m_alpha_x*m_alpha_x + sinPhi2*m_alpha_y*m_alpha_y ); 
	return alpha_i;
}

float MicrosurfaceSlopeBeckmann::P22(const float slope_x, const float slope_y) const
{
	const float value = 1.0f / (M_PI * m_alpha_x * m_alpha_y) * expf(-slope_x*slope_x/(m_alpha_x*m_alpha_x) - slope_y*slope_y/(m_alpha_y*m_alpha_y) );
	return value;
}

float MicrosurfaceSlopeBeckmann::Lambda(const vec3& wi) const
{
	if(wi.z > 0.9999f || wi.z < -0.9999f)
		return 0.0f;

	// a
	const float theta_i = acosf(wi.z);
	const float a = 1.0f/tanf(theta_i)/alpha_i(wi);

	// value
	const float value = 0.5f*((float)erf(a) - 1.0f) + INV_2_SQRT_M_PI / a * expf(-a*a);

	return value;
}

float MicrosurfaceSlopeBeckmann::projectedArea(const vec3& wi) const
{
	if(wi.z > 0.9999f)
		return 1.0f;
	if(wi.z < -0.9999f)
		return 0.0f;

	// a
	const float alphai = alpha_i(wi);
	const float theta_i = acosf(wi.z);
	const float a = 1.0f/tanf(theta_i)/alphai;

	// value
	const float value = 0.5f*((float)erf(a) + 1.0f)*wi.z + INV_2_SQRT_M_PI * alphai * sinf(theta_i) * expf(-a*a);

	return value;
}

vec2 MicrosurfaceSlopeBeckmann::sampleP22_11(const float theta_i, const float U, const float U_2) const
{
	vec2 slope;

	if(theta_i < 0.0001f)
	{
		const float r = sqrtf(-logf(U));
		const float phi = 6.28318530718f * U_2;
		slope.x = r * cosf(phi);
		slope.y = r * sinf(phi);
		return slope;
	}

	// constant
	const float sin_theta_i = sinf(theta_i);
	const float cos_theta_i = cosf(theta_i);

	// slope associated to theta_i
	const float slope_i = cos_theta_i/sin_theta_i;

	// projected area
	const float a = cos_theta_i/sin_theta_i;	
	const float projectedarea = 0.5f*((float)erf(a) + 1.0f)*cos_theta_i + INV_2_SQRT_M_PI * sin_theta_i * expf(-a*a);
	if(projectedarea < 0.0001f || projectedarea!=projectedarea)
		return vec2(0,0);
	// VNDF normalization factor
	const float c = 1.0f / projectedarea;

	// search 
	float erf_min = -0.9999f;
	float erf_max = std::max(erf_min, (float)erf(slope_i));
	float erf_current = 0.5f * (erf_min+erf_max);

	while(erf_max-erf_min > 0.00001f)
	{
		if (!(erf_current >= erf_min && erf_current <= erf_max))
			erf_current = 0.5f * (erf_min + erf_max);

		// evaluate slope
		const float slope = erfinv(erf_current);

		// CDF
		const float CDF = (slope>=slope_i) ? 1.0f : c * (INV_2_SQRT_M_PI*sin_theta_i*expf(-slope*slope) + cos_theta_i*(0.5f+0.5f*(float)erf(slope)));
		const float diff = CDF - U;

		// test estimate
		if( abs(diff) < 0.00001f )
			break;

		// update bounds
		if(diff > 0.0f)
		{
			if(erf_max == erf_current)
				break;
			erf_max = erf_current;
		}
		else
		{
			if(erf_min == erf_current)
				break;
			erf_min = erf_current;
		}

		// update estimate
		const float derivative = 0.5f*c*cos_theta_i - 0.5f*c*sin_theta_i * slope;
		erf_current -= diff/derivative;
	}

	slope.x = erfinv(std::min(erf_max, std::max(erf_min, erf_current)));
	slope.y = erfinv(2.0f*U_2-1.0f);
	return slope;
}

float MicrosurfaceSlopeGGX::P22(const float slope_x, const float slope_y) const
{
	const float tmp = 1.0f + slope_x*slope_x/(m_alpha_x*m_alpha_x) + slope_y*slope_y/(m_alpha_y*m_alpha_y);
	const float value = 1.0f / (M_PI * m_alpha_x * m_alpha_y) / (tmp * tmp);
	return value;
}

float MicrosurfaceSlopeGGX::Lambda(const vec3& wi) const
{
	if(wi.z > 0.9999f || wi.z < -0.9999f)
		return 0.0f;

	// a
	const float theta_i = acosf(wi.z);
	const float a = 1.0f/tanf(theta_i)/alpha_i(wi);

	// value
	const float value = 0.5f*(-1.0f + sign(a) * sqrtf(1 + 1/(a*a)));

	return value;
}

float MicrosurfaceSlopeGGX::projectedArea(const vec3& wi) const
{
	if(wi.z > 0.9999f)
		return 1.0f;
	if( wi.z < -0.9999f)
		return 0.0f;

	// a
	const float theta_i = acosf(wi.z);
	const float sin_theta_i = sinf(theta_i);

	const float alphai = alpha_i(wi);

	// value
	const float value = 0.5f * (wi.z + sqrtf(wi.z*wi.z + sin_theta_i*sin_theta_i*alphai*alphai));

	return value;
}

vec2 MicrosurfaceSlopeGGX::sampleP22_11(const float theta_i, const float U, const float U_2) const
{
	vec2 slope;

	if(theta_i < 0.0001f)
	{
		const float r = sqrtf(U/(1.0f-U));
		const float phi = 6.28318530718f * U_2;
		slope.x = r * cosf(phi);
		slope.y = r * sinf(phi);
		return slope;
	}

	// constant
	const float sin_theta_i = sinf(theta_i);
	const float cos_theta_i = cosf(theta_i);
	const float tan_theta_i = sin_theta_i/cos_theta_i;

	// slope associated to theta_i
	const float slope_i = cos_theta_i/sin_theta_i;

	// projected area
	const float projectedarea = 0.5f * (cos_theta_i + 1.0f);
	if(projectedarea < 0.0001f || projectedarea!=projectedarea)
		return vec2(0,0);
	// normalization coefficient
	const float c = 1.0f / projectedarea;

	const float A = 2.0f*U/cos_theta_i/c - 1.0f;
	const float B = tan_theta_i;
	const float tmp = 1.0f / (A*A-1.0f);

	const float D = sqrtf(std::max(0.0f, B*B*tmp*tmp - (A*A-B*B)*tmp));
	const float slope_x_1 = B*tmp - D;
	const float slope_x_2 = B*tmp + D;
	slope.x = (A < 0.0f || slope_x_2 > 1.0f/tan_theta_i) ? slope_x_1 : slope_x_2;

	float U2;
	float S;
	if(U_2 > 0.5f)
	{
	S = 1.0f;
	U2 = 2.0f*(U_2-0.5f);
	}
	else
	{
	S = -1.0f;
	U2 = 2.0f*(0.5f-U_2);
	}
	const float z = (U2*(U2*(U2*0.27385f-0.73369f)+0.46341f)) / (U2*(U2*(U2*0.093073f+0.309420f)-1.000000f)+0.597999f);
	slope.y = S * z * sqrtf(1.0f+slope.x*slope.x);

	return slope;
}


/************* MICROSURFACE *************/

float Microsurface::G_1(const vec3& wi) const
{
	if(wi.z > 0.9999f)
		return 1.0f;
	if(wi.z <= 0.0f)
		return 0.0f;

	// Lambda
	const float Lambda = m_microsurfaceslope->Lambda(wi);
	// value
	const float value = 1.0f / (1.0f + Lambda);
	return value;	
}

float Microsurface::G_1(const vec3& wi, const float h0) const
{
	if(wi.z > 0.9999f)
		return 1.0f;
	if(wi.z <= 0.0f)
		return 0.0f;

	// height CDF
	const float C1_h0 = m_microsurfaceheight->C1(h0);
	// Lambda
	const float Lambda = m_microsurfaceslope->Lambda(wi);
	// value
	const float value = powf(C1_h0, Lambda);
	return value;
}
float Microsurface::sampleHeight2(const vec3& wo, const float h0, const float U, float &throughput) const
{
	if(wo.z > 0.9999f)
		return FLT_MAX;
	if(wo.z < -0.9999f)
	{
		const float value = m_microsurfaceheight->invC1(U*m_microsurfaceheight->C1(h0));
		return value;
	}
	if(fabsf(wo.z) < 0.0001f)
		return h0;

	// probability of intersection
	const float G_1_ = G_1(wo, h0);
		
  throughput = 1.0-G_1_;
	return m_microsurfaceheight->invC1( 
			m_microsurfaceheight->C1(h0) / powf((1.0f-U*(1.0-G_1_)),1.0f/m_microsurfaceslope->Lambda(wo))
			);

}

float Microsurface::sampleHeight(const vec3& wo, const float h0, const float U) const
{
	if(wo.z > 0.9999f)
		return FLT_MAX;
	if(wo.z < -0.9999f)
	{
		const float value = m_microsurfaceheight->invC1(U*m_microsurfaceheight->C1(h0));
		return value;
	}
	if(fabsf(wo.z) < 0.0001f)
		return h0;

	// probability of intersection
	const float G_1_ = G_1(wo, h0);
		
	if (U >= 1.0f - G_1_) // leave the microsurface
		return FLT_MAX;

	const float h = m_microsurfaceheight->invC1( 
			m_microsurfaceheight->C1(h0) / powf((1.0f-U),1.0f/m_microsurfaceslope->Lambda(wo))
			);
	return h;
}

float Microsurface::sample(const vec3& wi, vec3& rayDir, int& scatteringOrder, float r0, float r1, float r2) const
{
	// init
	rayDir = -wi;
	float height_current = 1.0f + m_microsurfaceheight->invC1(0.999f);
	
	// random walk
	scatteringOrder = 0;	
  float throughput = 1.0f;
	while(true)
	{
		// next height
		float U = generateRandomNumber();
		height_current = sampleHeight(rayDir, height_current, U);

		// leave the microsurface?
		if( height_current == FLT_MAX )
			break;
		else
			scatteringOrder++;
    if(scatteringOrder > 5) return 0.0;

		// next direction
    vec3 newdir;
    if(scatteringOrder == 1)
      throughput *= samplePhaseFunction(-rayDir, newdir, r0, r1, r2);
    else
      throughput *= samplePhaseFunction(-rayDir, newdir, generateRandomNumber(), generateRandomNumber(), generateRandomNumber());
    rayDir = newdir;

		// if NaN (should not happen, just in case ;-)
		if(!(height_current == height_current)) return 0.0f;
	}

	return throughput;
}

float Microsurface::eval(const vec3& wi, const vec3& wo, const int scatteringOrder) const
{
  if(fabsf(wo[2]) < 1e-5f) return 0.0f;
	// init
	vec3 rayDir = -wi;
	float height_current = 1.0f + m_microsurfaceheight->invC1(0.999f);

	float sum = 0;
  float throughput = 1.0f;
	
	// random walk
	int current_scatteringOrder = 0;	
	while(scatteringOrder==0 || current_scatteringOrder <= scatteringOrder)
	{
		// next height
		float U = generateRandomNumber();
    // float thr = 1.0f;
		// height_current = sampleHeight2(rayDir, height_current, U, thr);
    // throughput *= thr;
		height_current = sampleHeight(rayDir, height_current, U);

		// leave the microsurface?
		if( height_current == FLT_MAX )
			break;
		else
			current_scatteringOrder++;
    if(current_scatteringOrder > 5) return sum;

		// next event estimation
		float phasefunction = evalPhaseFunction(-rayDir, wo);
		float shadowing = G_1(wo, height_current);
		float I = phasefunction * shadowing;

		// XXX if ( IsFiniteNumber(I) && (scatteringOrder==0 || current_scatteringOrder==scatteringOrder) )
		if ( IsFiniteNumber(I))// && (current_scatteringOrder>1))
			sum += throughput * I;
		
		// next direction
    vec3 newdir;
    throughput *= samplePhaseFunction(-rayDir, newdir, generateRandomNumber(), generateRandomNumber(), generateRandomNumber());
    rayDir = newdir;
			
		// if NaN (should not happen, just in case ;-)
		if(height_current != height_current) return 0.0f;
	}

	return sum/fabsf(wo[2]);
}
namespace{
float CondFresnel(const float n1, const std::complex<float> n2, float cosr)
{
  cosr = std::max(0.0f, cosr);
  const complex<float> cost = std::sqrt(1.0f - (n1/n2)*(n1/n2) * (1.0f - cosr*cosr));
  // fresnel for unpolarized light:
  const float Rs = std::abs((n1*cosr - n2*cost)/(n1*cosr + n2*cost));
  const float Rp = std::abs((n1*cost - n2*cosr)/(n1*cost + n2*cosr));
  return std::min(1.0f, (Rs*Rs + Rp*Rp)*.5f);
}
}

float MicrosurfaceConductor::evalPhaseFunction(const vec3& wi, const vec3& wo) const
{
	// half vector 
	const vec3 wh = normalize(wi+wo);
	if(wh.z < 0.0f)
		return 0.0f;

  const float dot_wi_wh = dot(wi, wh);
  if(std::abs(dot_wi_wh) < 1e-5f) return 0.0f;
	
  const float F = CondFresnel(1.0, m_ior, dot_wi_wh);
	// value
	const float value = 0.25f * F * m_microsurfaceslope->D_wi(wi, wh) / dot_wi_wh;
	return value;
}

float MicrosurfaceConductor::samplePhaseFunction(const vec3& wi, vec3 &wo, float U0, float U1, float U2) const
{
	vec3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

	// reflect
	wo = -wi + 2.0f * wm * dot(wi, wm);

	return CondFresnel(1.0, m_ior, dot(wi, wm));
}

float MicrosurfaceConductor::evalSingleScattering(const vec3& wi, const vec3& wo) const
{
	// half-vector
	const vec3 wh = normalize(wi+wo);
	const float D = m_microsurfaceslope->D(wh);

	// masking-shadowing 
	const float G2 = 1.0f / (1.0f + m_microsurfaceslope->Lambda(wi) + m_microsurfaceslope->Lambda(wo));

  // fresnel
  const float F = CondFresnel(1.0, m_ior, dot(wi, wh));

	// BRDF * cos
	const float value = F * D * G2 / (4.0f * wi.z);

	return value;
}

static float fresnelDielectricExt(float cosThetaI_, float &cosThetaT_, float eta) 
{

	/* Using Snell's law, calculate the squared sine of the
	   angle between the normal and the transmitted ray */
	float scale = (cosThetaI_ > 0) ? 1/eta : eta,
	      cosThetaTSqr = 1 - (1-cosThetaI_*cosThetaI_) * (scale*scale);

	/* Check for total internal reflection */
	if (cosThetaTSqr <= 0.0f) {
		cosThetaT_ = 0.0f;
		return 1.0f;
	}

	/* Find the absolute cosines of the incident/transmitted rays */
	float cosThetaI = std::abs(cosThetaI_);
	float cosThetaT = std::sqrt(cosThetaTSqr);

	float Rs = (cosThetaI - eta * cosThetaT)
			 / (cosThetaI + eta * cosThetaT);
	float Rp = (eta * cosThetaI - cosThetaT)
			 / (eta * cosThetaI + cosThetaT);

	cosThetaT_ = (cosThetaI_ > 0) ? -cosThetaT : cosThetaT;

	/* No polarization -- return the unpolarized reflectance */
	return 0.5f * (Rs * Rs + Rp * Rp);
}

static vec3 refract(const vec3 &wi, const vec3 &n, float eta, float cosThetaT)  
{
	if (cosThetaT < 0)
		eta = 1.0f / eta;
	return n * (dot(wi, n) * eta + cosThetaT) - wi * eta;
}


float MicrosurfaceDielectric::Fresnel(const vec3& wi, const vec3& wm, const float eta) const
{
	float cosThetaT_;
	return fresnelDielectricExt(dot(wi, wm), cosThetaT_, eta);
}

float MicrosurfaceDielectric::evalPhaseFunction(const vec3& wi, const vec3& wo) const
{
	return evalPhaseFunction(wi, wo, true, true) + evalPhaseFunction(wi, wo, true, false);
}

float MicrosurfaceDielectric::evalPhaseFunction(const vec3& wi, const vec3& wo, const bool wi_outside, const bool wo_outside) const
{
	const float eta = wi_outside ? m_eta : 1.0f / m_eta;

	if( wi_outside == wo_outside ) // reflection
	{
    // fprintf(stderr, "outside %d %d wi wo %g %g\n", wi_outside, wo_outside, wi.z, wo.z);
		// half vector 
		const vec3 wh = normalize(wi+wo);
		// value
		const float value = (wi_outside) ?
								(0.25f * m_microsurfaceslope->D_wi(wi, wh) / dot(wi, wh) * Fresnel(wi, wh, eta)) :
								(0.25f * m_microsurfaceslope->D_wi(-wi, -wh) / dot(-wi, -wh) * Fresnel(-wi, -wh, eta)) ;
    // if(value > 0)
    // {
    // if(wi_outside) fprintf(stderr, "eval %g D_wi %g %g\n", value, wi.z, wh.z);
    // else fprintf(stderr, "eval %g D_wi %g %g\n", value, -wi.z, -wh.z);
    // }
		return value;
	}
	else // transmission
	{
		vec3 wh = -normalize(wi+wo*eta);
		wh *= (wi_outside) ? (sign(wh.z)) : (-sign(wh.z));

		if(dot(wh, wi) < 0)
			return 0;

		const float value = (wi_outside) ? 
			(1.0f-Fresnel(wi, wh, eta)) * m_microsurfaceslope->D_wi(wi, wh) * eta*eta * std::max(0.0f, -dot(wo, wh)) / powf(dot(wi, wh)+eta*dot(wo,wh), 2.0f) :
			(1.0f-Fresnel(-wi, -wh, eta)) * m_microsurfaceslope->D_wi(-wi, -wh) * eta*eta * std::max(0.0f, -dot(-wo, -wh)) / powf(dot(-wi, -wh)+eta*dot(-wo,-wh), 2.0f);
		return value;	
	}
}

float MicrosurfaceDielectric::samplePhaseFunction(const vec3& wi, vec3& wo, const bool wi_outside, bool& wo_outside,
    float U0, float U1, float U2) const
{
	const float eta = wi_outside ? m_eta : 1.0f / m_eta;

	vec3 wm = wi_outside ? (m_microsurfaceslope->sampleD_wi(wi, U1, U2)) :
						   (-m_microsurfaceslope->sampleD_wi(-wi, U1, U2)) ;

	float cosThetaT;
	const float F = fresnelDielectricExt(dot(wi, wm), cosThetaT, eta);

	if( U0 < F )
	{
		wo = -wi + 2.0f * wm * dot(wi, wm); // reflect
		return 1.0f;
	}
	else
	{
		wo_outside = !wi_outside;
		wo = normalize(refract(wi, wm, eta, cosThetaT));
    return 1.0f;
	}
}

float MicrosurfaceDielectric::evalSingleScattering(const vec3& wi, const vec3& wo) const
{
	bool wi_outside = true;
	bool wo_outside = wo.z > 0;

	const float eta = m_eta;

	if(wo_outside) // reflection
	{
		// D
		const vec3 wh = normalize(vec3(wi+wo));
		const float D = m_microsurfaceslope->D(wh);
		
		// masking shadowing
		const float Lambda_i = m_microsurfaceslope->Lambda(wi);
		const float Lambda_o = m_microsurfaceslope->Lambda(wo);
		const float G2 = 1.0f / (1.0f + Lambda_i + Lambda_o);

		// BRDF
		const float value = Fresnel(wi, wh, eta) * D * G2 / (4.0f * wi.z);
		return value;
	}
	else // refraction
	{
		// D
		vec3 wh = -normalize(wi+wo*eta);
		if(eta<1.0f)
			wh = -wh;
		const float D = m_microsurfaceslope->D(wh);

		// G2
		const float Lambda_i = m_microsurfaceslope->Lambda(wi);
		const float Lambda_o = m_microsurfaceslope->Lambda(-wo);
		const float G2 = (float) beta(1.0f+Lambda_i, 1.0f+Lambda_o);

		// BSDF
		const float value = std::max(0.0f, dot(wi, wh)) * std::max(0.0f, -dot(wo, wh)) / wi.z * eta*eta * G2 * D / powf(dot(wi, wh)+eta*dot(wo,wh), 2.0f) * (1.0f-Fresnel(wi, wh, eta));
		return value;
	}
}

float MicrosurfaceDielectric::eval(const vec3& wi, const vec3& wo, const int scatteringOrder) const
{
  if(fabsf(wo[2]) < 1e-5f) return 0.0f;
	// init
	vec3 rayDir = -wi;
	float height_current = 1.0f + m_microsurfaceheight->invC1(0.999f);
	bool outside = true;

	float sum = 0.0f;
	
	// random walk
	int current_scatteringOrder = 0;	
  float throughput = 1.0f;
	while(scatteringOrder==0 || current_scatteringOrder <= scatteringOrder)
	{
		// next height
		float U = generateRandomNumber();		
    // float thr = 1.0f;
		// height_current = (outside) ? sampleHeight2(rayDir, height_current, U, thr) : -sampleHeight2(-rayDir, -height_current, U, thr);
    // throughput *= thr;
		height_current = (outside) ? sampleHeight(rayDir, height_current, U) : -sampleHeight(-rayDir, -height_current, U);

		// leave the microsurface?
		if( height_current == FLT_MAX || height_current == -FLT_MAX)
			break;
		else
			current_scatteringOrder++;
    if(current_scatteringOrder > 5) return sum;

		// next event estimation
		float phasefunction = evalPhaseFunction(-rayDir, wo, outside, (wo.z>0) );
		float shadowing = (wo.z>0) ? G_1(wo, height_current) : G_1(-wo, -height_current);
		float I = phasefunction * shadowing;

		// if ( IsFiniteNumber(I) && (scatteringOrder==0 || current_scatteringOrder==scatteringOrder) )
		// if ( IsFiniteNumber(I) && (current_scatteringOrder>1) ) // XXX kill direct
		if ( IsFiniteNumber(I))
			sum += throughput * I;
		
		// next direction
    vec3 newdir;
		throughput *= samplePhaseFunction(-rayDir, newdir, outside, outside, generateRandomNumber(), generateRandomNumber(), generateRandomNumber());
    rayDir = newdir;

		// if NaN (should not happen, just in case ;-)
		if(height_current != height_current) return 0.0f;
	}

	return sum / fabsf(wo[2]);
}

#if 0
float Microsurface::eval_precached(const vec3& wi, const vec3& wo, const vertex_multiscatter_t *s, const int bounces) const
{
  return 0.;
}
float MicrosurfaceDielectric::eval_precached(const vec3& wi, const vec3& wo, const vertex_multiscatter_t *s, const int bounces) const
{
  // TODO: first bounce separately
  float sum = 0.0f;
  for(int i=0;i<bounces;i++)
  {
		// next event estimation
		float phasefunction = evalPhaseFunction(-vec3(s[i].wi[0], s[i].wi[1], s[i].wi[2]), wo, !s[i].inside, (wo.z>0) );
		float shadowing = (wo.z>0) ? G_1(wo, s[i].height) : G_1(-wo, -s[i].height);
		float I = s[i].throughput * phasefunction * shadowing;

		// if ( IsFiniteNumber(I) && (scatteringOrder==0 || current_scatteringOrder==scatteringOrder) )
		if ( IsFiniteNumber(I))
			sum += I;
  }
  return sum;
}
#endif

float Microsurface::pdf(const vec3& wi, const vec3& wo) const
{
  if(std::abs(wo[2]) < 1e-5f) return 0.0f;
  return 0.5f/M_PI+.5f*evalPhaseFunction(wi, wo)/std::abs(wo[2]);
}
float MicrosurfaceDielectric::pdf(const vec3& wi, const vec3& wo) const
{
  if(std::abs(wo[2]) < 1e-5f) return 0.0f;
  return 0.5f/M_PI+ .5f*evalPhaseFunction(wi, wo, true, (wo.z>0))/std::abs(wo[2]);
}

#if 0
int Microsurface::precache(const vec3& wi, vertex_multiscatter_t *s) const
{
  return 0;
}
int MicrosurfaceDielectric::precache(const vec3& wi, vertex_multiscatter_t *s) const
{
	// init
	vec3 rayDir = -wi;
	float height_current = 1.0f + m_microsurfaceheight->invC1(0.999f);
	bool outside = true;
	
	// random walk
	int scatteringOrder = 0;	
  int caches = 0;
  float throughput = 1.0f;
	while(true)
	{
		// next height
		float U = generateRandomNumber();

		height_current = (outside) ? sampleHeight(rayDir, height_current, U) : -sampleHeight(-rayDir, -height_current, U);

		// leave the microsurface?
		if( height_current == FLT_MAX || height_current == -FLT_MAX)
      return caches;
		else
			scatteringOrder++;
    if(scatteringOrder > 10) return caches;

    // if(scatteringOrder > 1) // first bounce separately
    {
      int i = caches++;
      if(i < MICRO_MAX_BOUNCES)
      {
        s[i].height = height_current;
        s[i].inside = !outside;
        s[i].wi[0] = rayDir.x;
        s[i].wi[1] = rayDir.y;
        s[i].wi[2] = rayDir.z;
        s[i].throughput = throughput;
      }
    }

		// next direction
    vec3 newdir;
		throughput *= samplePhaseFunction(-rayDir, newdir, outside, outside, generateRandomNumber(), generateRandomNumber(), generateRandomNumber());
    rayDir = newdir;

		// if NaN (should not happen, just in case ;-)
		if(height_current != height_current) return caches;
	}

  // never reached
	return 0;
}
#endif


float MicrosurfaceDielectric::sample(const vec3& wi, vec3 &wo, int& scatteringOrder, float r0, float r1, float r2) const
{
	// init
	wo = -wi;
	float height_current = 1.0f + m_microsurfaceheight->invC1(0.999f);
	bool outside = true;
	
	// random walk
	scatteringOrder = 0;	
  float throughput = 1.0f;
	while(true)
	{
		// next height
		float U = generateRandomNumber();

		height_current = (outside) ? sampleHeight(wo, height_current, U) : -sampleHeight(-wo, -height_current, U);

		// leave the microsurface?
		if( height_current == FLT_MAX || height_current == -FLT_MAX)
    {
      // XXX kill direct first bounce:
      // if(scatteringOrder <= 1) return 0.0f;
			break;
    }
		else
			scatteringOrder++;
    if(scatteringOrder > 5) return 0;

		// next direction
    vec3 newdir;
    if(scatteringOrder == 1)
      throughput *= samplePhaseFunction(-wo, newdir, outside, outside, r0, r1, r2);
    else
      throughput *= samplePhaseFunction(-wo, newdir, outside, outside, generateRandomNumber(), generateRandomNumber(), generateRandomNumber());
    wo = newdir;

		// if NaN (should not happen, just in case ;-)
		if(height_current != height_current) return 0.0f;
	}

	return throughput;
}

float MicrosurfaceDiffuse::evalPhaseFunction(const vec3& wi, const vec3& wo) const
{
	const float U1 = generateRandomNumber();
	const float U2 = generateRandomNumber();
	vec3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);
	
	// value
	const float value = 1.0f/M_PI * std::max(0.0f, dot(wo, wm));
	return value * m_reflectivity;
}

 // build orthonormal basis (Building an Orthonormal Basis from a 3D Unit Vector Without Normalization, [Frisvad2012])
 void buildOrthonormalBasis(vec3& omega_1, vec3& omega_2, const vec3& omega_3)
{
	if(omega_3.z < -0.9999999f) 
	{
	   omega_1 = vec3 ( 0.0f , -1.0f , 0.0f );
	   omega_2 = vec3 ( -1.0f , 0.0f , 0.0f );
	} else {
	   const float a = 1.0f /(1.0f + omega_3.z );
	   const float b = -omega_3.x*omega_3 .y*a ;
	   omega_1 = vec3 (1.0f - omega_3.x*omega_3. x*a , b , -omega_3.x );
	   omega_2 = vec3 (b , 1.0f - omega_3.y*omega_3.y*a , -omega_3.y );
	}
}


float MicrosurfaceDiffuse::samplePhaseFunction(const vec3& wi, vec3& out, float U1, float U2, float U3) const
{
	const float U4 = generateRandomNumber();

	vec3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

	// sample diffuse reflection
	vec3 w1, w2;
	buildOrthonormalBasis(w1, w2, wm);

	float r1 = 2.0f*U3 - 1.0f;
	float r2 = 2.0f*U4 - 1.0f;

	// concentric map code from
	// http://psgraphics.blogspot.ch/2011/01/improved-code-for-concentric-map.html
	float phi, r;
	if (r1 == 0 && r2 == 0) {
		r = phi = 0;
	} else if (r1*r1 > r2*r2) {
		r = r1;
		phi = (M_PI/4.0f) * (r2/r1);
	} else {
		r = r2;
		phi = (M_PI/2.0f) - (r1/r2) * (M_PI/4.0f);
	}
	float x = r*cosf(phi);
	float y = r*sinf(phi);
	float z = sqrtf(1.0f - x*x - y*y);
	out = x*w1 + y*w2 + z*wm;
  return m_reflectivity;
}

// stochastic evaluation  
// Heitz and Dupuy 2015
// Implementing a Simple Anisotropic Rough Diffuse Material with Stochastic Evaluation
float MicrosurfaceDiffuse::evalSingleScattering(const vec3& wi, const vec3& wo) const
{
	// sample visible microfacet
	const float U1 = generateRandomNumber();
	const float U2 = generateRandomNumber();
	const vec3 wm = m_microsurfaceslope->sampleD_wi(wi, U1, U2);

	// shadowing given masking
	const float Lambda_i = m_microsurfaceslope->Lambda(wi);
	const float Lambda_o = m_microsurfaceslope->Lambda(wo);
	float G2_given_G1 = (1.0f + Lambda_i) / (1.0f + Lambda_i + Lambda_o);

	// evaluate diffuse and shadowing given masking
	const float value = 1.0f / (float)M_PI * std::max(0.0f, dot(wm, wo)) * G2_given_G1;

	return m_reflectivity * value;
}

