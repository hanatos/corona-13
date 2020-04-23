#include <math.h>
#include <fftw3.h>
#ifndef M_PI
# define M_PI           3.14159265358979323846  /* pi */
#endif

#define NUM_THREADS 12 //set number of threads to number of available processors

// based on pravin's screened 2d poisson solver http://grail.cs.washington.edu/projects/screenedPoissonEq/screenedPoissonEq_files/fft_solver.cpp

/*
 * The following 100 line function is our entire FFT based solver for the screened Poisson equation.
 * This function won't compile without our Image class (not included but you could replace it with your own) and the fftw library (http://www.fftw.org/). 
 */


static inline void fourier_solve(
    float *img,  // will be overwritten
    const float *grad_x,
    const float *grad_y,
    const int width,
    const int height,
    const float alpha)
{
  fftwf_init_threads();

  float* fftBuff   = (float*) fftwf_malloc(sizeof(*fftBuff)   * width * height);

  //compute two 1D lookup tables for computing the DCT of a 2D Laplacian on the fly
  float* ftLapY = (float*) fftwf_malloc(sizeof(*ftLapY) * height);
  float* ftLapX = (float*) fftwf_malloc(sizeof(*ftLapX) * width);

  for(int x = 0; x < width; x++)
  {
    ftLapX[x] = 2.0f * cosf(M_PI * x / (width - 1));
  }
  for(int y = 0; y < height; y++)
  {
    ftLapY[y] = -4.0f + (2.0f * cosf(M_PI * y / (height - 1)));
  }

  //Create a DCT-I plan for, which is its own inverse.
  fftwf_plan_with_nthreads(NUM_THREADS);
  fftwf_plan fftPlan;
  // fftPlan = fftwf_plan_r2r_2d(height, width, fftBuff, fftBuff, 
  fftPlan = fftwf_plan_r2r_2d(height, width, 0, 0,
      FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE); //use FFTW_PATIENT when plan can be reused

  for(int c = 0; c < 3; c++)
  {
    printf("solving channel %i\n", c); 

    int nodeAddr        = 0;
    int pixelAddr       = c;
    int rightPixelAddr  = 3 + c;
    int topPixelAddr    = (width * 3) + c;

    float dcSum = 0.0f;

    // compute h_hat from u, gx, gy (see equation 48 in the paper), as well as the DC term of u's DCT.
    for(int y = 0; y < height; y++)
      for(int x = 0; x < width;  x++, 
          nodeAddr++, pixelAddr += 3, rightPixelAddr += 3, topPixelAddr += 3)
      {
        // Compute DC term of u's DCT without computing the whole DCT.
        float dcMult = 1.0f;
        if((x > 0) && (x < width  - 1))
          dcMult *= 2.0f;
        if((y > 0) && (y < height - 1))
          dcMult *= 2.0f;
        dcSum += dcMult * img[pixelAddr];

        fftBuff[nodeAddr] = alpha * img[pixelAddr];

        // Subtract g^x_x and g^y_y, with boundary factor of -2.0 to account for boundary reflections implicit in the DCT
        if((x > 0) && (x < width - 1))
          fftBuff[nodeAddr] -= grad_x[rightPixelAddr] - grad_x[pixelAddr];
        else
          fftBuff[nodeAddr] -= -2.0f * grad_x[pixelAddr];

        if((y > 0) && (y < height - 1))
          fftBuff[nodeAddr] -= grad_y[topPixelAddr] - grad_y[pixelAddr];
        else
          fftBuff[nodeAddr] -= -2.0f * grad_y[pixelAddr];
      }

    //transform h_hat to H_hat by taking the DCT of h_hat
    fftwf_execute(fftPlan);

    //compute F_hat using H_hat (see equation 29 in the paper)
    nodeAddr = 0;
    for(int y = 0; y < height; y++)
      for(int x = 0; x < width;  x++, nodeAddr++)
      {
        float ftLapResponse = ftLapY[y] + ftLapX[x]; 
        fftBuff[nodeAddr] /= (alpha - ftLapResponse);
      }

    /*
     * Set the DC term of the solution to the value computed above (i.e., the DC term of imgData). 
     * When dataCost = 0 (i.e., there is no data image and the problem becomes pure gradient field integration)
     * then the DC term  of the solution is undefined. So if you want to control the DC of the solution 
     * when dataCost = 0 then before calling fourierSolve() set every pixel in 'imgData' to the average value 
     * you would like the pixels in the solution to have. 
     */
    fftBuff[0] = dcSum;

    //transform F_hat to f_hat by taking the inverse DCT of F_hat
    fftwf_execute(fftPlan);
    float fftDenom = 4.0f * (width - 1) * (height - 1);
    pixelAddr = c;
    for(int k = 0; k < width*height; k++, pixelAddr += 3)
    {
      img[pixelAddr] = fftBuff[k] / fftDenom;
    }
  }

  fftwf_free(fftBuff);
  fftwf_free(ftLapX);
  fftwf_free(ftLapY);
  fftwf_destroy_plan(fftPlan);
  fftwf_cleanup_threads();

  printf("done.\n");
}

