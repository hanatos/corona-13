// ==========================================================================
// $Id$ 
// Interface class for multispectral BRDF. Reads from file and performs
// "quintilinear" interpolation
// ==========================================================================
// (C)opyright:
//
//   Max-Planck-Institut fuer Informatik
//   Campus E1.4, 66123 Saarbruecken
//   http://www.mpi-sb.mpg.de
// 
// Creator: hullin (Matthias Hullin)
// Email:   hullin@mpi-sb.mpg.de
// ==========================================================================
// $Log$
// ==========================================================================
#ifndef MULTISPECBRDF_H
#define MULTISPECBRDF_H

#define sqr(x) ((x)*(x))
#define PI 3.14159265358979323846

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <assert.h>

/*! \file multispecBRDF.hh 
    \brief Interface class for multispectral
    BRDF. Reads from file and performs "quintilinear" interpolation
 */


  /** \class multispecBRDF multispecBRDF.hh Interface class for
      multispectral BRDF. Reads from file and performs "quintilinear"
      interpolation */

using namespace std;


  
  class multispecBRDF {

  public:
#ifdef REDUCE
   int reduced;
#endif

    /** default destructor */
    ~multispecBRDF() { delete[](data); }

    /** default constructor; please avoid */
    multispecBRDF() {};

    /** init from buffer constructor. */
    multispecBRDF(char *buf);

    /** "load file" constructor */
    multispecBRDF(const string &filename) ;



    /** compute linear index, given a multidim set of bin indices */
    inline size_t index(int theta_in_bin, 
			int theta_out_bin, 
			int phi_rel_bin, 
			int lambda_in_bin, 
			int lambda_out_bin) {
      return (  (size_t)lambda_out_bin * dir[4]
	      + (size_t)lambda_in_bin * dir[3]
	      + (size_t)phi_rel_bin  * dir[2]
	      + (size_t)theta_out_bin * dir[1]
	      + (size_t)theta_in_bin * dir[0] );
    };


#ifdef REDUCE
    void reduce(){
	static int _reduced = 0;
	    // reduced=true;
	    reduced = ++_reduced;
	    dir[4] <<= 1;
	    dir[3] <<= 1;
	    _illumStep <<= 1;
	    _camStep <<= 1;
	    _illumBins >>= 1;
	    _camBins >>=1;

	    for (int i=0; i<32; ++i) {
		    size_t b0 = (i & 16) >> 4;
		    size_t b1 = (i & 8) >> 3;
		    size_t b2 = (i & 4) >> 2;
		    size_t b3 = (i & 2) >> 1;
		    size_t b4 = (i & 1);
		    dircube[i] = b0 * dir[0] + b1 * dir[1] + b2 * dir[2] + b3 * dir[3] + b4 * dir[4];
	    }
    }; 
#endif

    void smoothSpecularHighlight(int  boxHalfWidth = 2, int phiTolerance = 3) {


	    float *temp = new float[ _thetaBins * _thetaBins * phiTolerance ];


	    for (int lo = 0; lo < _camBins; ++lo) {
		    int li = lo + _illumBins - _camBins;
		    // cout <<"i"<<li<<"o"<<lo<<endl;
		    for (int p = _phiBins - phiTolerance; p < _phiBins; ++p) {
			    // cout <<p<<" "<<(_phiBins - 1 - p)<<endl;
			    for (int ti = 0; ti < _thetaBins; ++ti)
				    for (int to = 0; to < _thetaBins; ++to) {
					    float sum = 0;
					    int entries = 0;

					    for (int i = -boxHalfWidth; i <= boxHalfWidth; ++i)
						    if (ti + i >= 0 && ti + i < _thetaBins
								    && to + i >= 0 && to + i < _thetaBins) {
							    sum += data[index( ti + i, to + i, p, li, lo)];
							    entries += 1;
						    }

					    sum /= entries;
					    temp[ti + to *_thetaBins + (_phiBins - 1 - p) * _thetaBins * _thetaBins] = sum;
				    }
		    }

		    for (int p = _phiBins -3; p<_phiBins; ++p)
			    for (int ti=0; ti<_thetaBins; ++ti)
				    for (int to=0; to<_thetaBins; ++to)
					    data[index( ti, to, p, li, lo)] = temp[ti + to *_thetaBins + (_phiBins - 1 - p) * _thetaBins * _thetaBins];

	    }

	    delete[] temp;
    } 


    /** interpolated value readout. Believe it or not, a 5-fold loop
	is the absolute maximum of what you want to unroll :-) */

    inline float interpolate(float theta_in_bin, 
			     float theta_out_bin, 
			     float phi_rel_bin, 
			     float lambda_in_bin, 
			     float lambda_out_bin) {
      
      float idx[5] = {theta_in_bin, theta_out_bin, phi_rel_bin, lambda_in_bin, lambda_out_bin};
      
      int idx_int[5] = {(int)floor(idx[0]), (int)floor(idx[1]), (int)floor(idx[2]), (int)floor(idx[3]), (int)floor(idx[4])}; 
      
      float idx_frac[5] = {idx[0] - idx_int[0], 
			   idx[1] - idx_int[1], 
			   idx[2] - idx_int[2], 
			   idx[3] - idx_int[3], 
			   idx[4] - idx_int[4]};

      // avoid range violation on data readout
      // if (idx_int[0] == _thetaBins-1) {idx_int[0] = _thetaBins-2; idx_frac[0] = 1; }
      // else // this "else" skips quintuple check in most cases
      // {
      //   if (idx_int[1] == _thetaBins-1) {idx_int[1] = _thetaBins-2; idx_frac[1] = 1; }
      //   if (idx_int[2] == _phiBins-1) {idx_int[2] = _phiBins-2; idx_frac[2] = 1; }
      //   if (idx_int[3] == _illumBins-1) {idx_int[3] = _illumBins-2; idx_frac[3] = 1; }
      //   if (idx_int[4] == _camBins-1) {idx_int[4] = _camBins-2; idx_frac[4] = 1; }
      // }
      if (idx_int[0] > _thetaBins - 2) {idx_int[0] = _thetaBins - 2; idx_frac[0] = 1; }
      
      if (idx_int[1] > _thetaBins - 2) {idx_int[1] = _thetaBins - 2; idx_frac[1] = 1; }
      
      if (idx_int[2] > _phiBins - 2) {idx_int[2] = _phiBins - 2; idx_frac[2] = 1; }

      size_t offset = index(idx_int[0], idx_int[1], idx_int[2], idx_int[3], idx_int[4]);

      float v[32]; // holds a data hypercube
      for (int i=0; i<32; ++i) {
	v[i] = data[offset + dircube[i]];
      }
      
      float b = idx_frac[0];
      float mb = 1-idx_frac[0];
	
      v[0] = v[0] * mb + v[16] * b;
      v[1] = v[1] * mb + v[17] * b;
      v[2] = v[2] * mb + v[18] * b;
      v[3] = v[3] * mb + v[19] * b;
      v[4] = v[4] * mb + v[20] * b;
      v[5] = v[5] * mb + v[21] * b;
      v[6] = v[6] * mb + v[22] * b;
      v[7] = v[7] * mb + v[23] * b;
      v[8] = v[8] * mb + v[24] * b;
      v[9] = v[9] * mb + v[25] * b;
      v[10] = v[10] * mb + v[26] * b;
      v[11] = v[11] * mb + v[27] * b;
      v[12] = v[12] * mb + v[28] * b;
      v[13] = v[13] * mb + v[29] * b;
      v[14] = v[14] * mb + v[30] * b;
      v[15] = v[15] * mb + v[31] * b;

      b = idx_frac[1];
      mb = 1-idx_frac[1];
      v[0] = v[0] * mb + v[8] * b;
      v[1] = v[1] * mb + v[9] * b;
      v[2] = v[2] * mb + v[10] * b;
      v[3] = v[3] * mb + v[11] * b;
      v[4] = v[4] * mb + v[12] * b;
      v[5] = v[5] * mb + v[13] * b;
      v[6] = v[6] * mb + v[14] * b;
      v[7] = v[7] * mb + v[15] * b;

      b = idx_frac[2];
      mb = 1-idx_frac[2];
      v[0] = v[0] * mb + v[4] * b;
      v[1] = v[1] * mb + v[5] * b;
      v[2] = v[2] * mb + v[6] * b;
      v[3] = v[3] * mb + v[7] * b;

      b = idx_frac[3];
      mb = 1-idx_frac[3];
      v[0] = v[0] * mb + v[2] * b;
      v[1] = v[1] * mb + v[3] * b;

      return v[0] * (1-idx_frac[4]) + v[1] * idx_frac[4];
     
    };

    /* interpolated readout accepting arguments in real-world units */
    inline float readout(float theta_in,  // range [0..Pi/2]
			float theta_out,    // range [0..Pi/2]
			float phi_rel,      // range [0..Pi]
			float lambda_in,   // all wavelengths in nm
			float lambda_out) {
#if 0 // stupidly stretch theta to < 80 deg.
      theta_in *= 0.8888f;
      theta_out *= 0.8888f;
#endif
      
      float p = _phiBins*phi_rel/PI-0.5;
      float ti = _thetaBins*theta_in*2/PI-0.5;
      float to = _thetaBins * theta_out * 2 / PI - 0.5;
      

      if (ti<0.f) ti=0.f; 
      if (to<0.f) to=0.f; 
      if (p<0.f) p=0.f; 
      // float to = _thetaBins*theta_out*2/PI-0.5;
      // if (!(ti>=0.f)) ti=0.f; else if (ti>_thetaBins-1.f) ti = _thetaBins - 1.00001f;
      // if (!(to>=0.f)) to=0.f; else if (to>_thetaBins-1.f) to = _thetaBins - 1.00001f;
      // if (!(p >=0.f)) p=0.f; else if (p>_phiBins-1.f) p = _phiBins - 1.00001f;

#ifdef DOWNSAMPLE
      // downsample lambda in/out in bins
      lambda_in  = 400.0f + (.5f+(int)((lambda_in  - 400.0f)*DOWNSAMPLE/300.0f))*300.0f/DOWNSAMPLE;
      lambda_out = 400.0f + (.5f+(int)((lambda_out - 400.0f)*DOWNSAMPLE/300.0f))*300.0f/DOWNSAMPLE;
#endif
#if 0
      int li = (lambda_in-_illumFrom)/_illumStep;
      int lo = (lambda_out-_camFrom)/_camStep;
      if (!(li>=0)) li=0; else if (li>_illumBins-1) li = _illumBins - 1;
      if (!(lo>=0)) lo=0; else if (lo>_camBins-1) lo = _camBins - 1;
#else
      float li = (lambda_in-_illumFrom)/_illumStep;
      float lo = (lambda_out-_camFrom)/_camStep;
      if (!(li>=0.f)) li=0.f; // else if (li>_illumBins-1.f) li = _illumBins - 1.00001f;
      if (!(lo>=0.f)) lo=0.f; // else if (lo>_camBins-1.f) lo = _camBins - 1.00001f;
#endif
#ifdef DIAG
      if((int)li != (int)lo) return .0f;
#endif
#ifdef REDUCE
      return (1<<reduced)*(1<<reduced)*interpolate(ti, to, p, li, lo)/_illumStep;
#else
      // return interpolate(ti, to, p, li, lo)/_illumStep;
      return interpolate(ti, to, p, li, lo); // the division by _illumStep is necessary for dimensions, see paper
#endif
    }

    /* give reference to a bin */
    inline float & operator ()(int theta_in_bin, 
			int theta_out_bin, 
			int phi_rel_bin, 
			int lambda_in_bin, 
			int lambda_out_bin) {
      return data[index( theta_in_bin, theta_out_bin, phi_rel_bin, lambda_in_bin, lambda_out_bin)];
    };

  protected:
    float* data;
    int _thetaBins, _phiBins;
    int _illumBins, _camBins;
    int _illumFrom, _illumStep, _illumTo; // wavelength ranges (light/camera) of measurement
    int _camFrom, _camStep, _camTo;

    size_t dir[5];
    size_t dircube[32];

  };





#endif /* MULTISPECBRDF_H */







