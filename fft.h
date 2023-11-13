#ifndef _FFT_H_
#define _FFT_H_

/**********************************************************************

  FFT.h

  Dominic Mazzoni

  September 2000

  This file contains a few FFT routines, including a real-FFT
  routine that is almost twice as fast as a normal complex FFT,
  and a power spectrum routine which is more convenient when
  you know you don't care about phase information.  It now also
  contains a few basic windowing functions.

  Some of this code was based on a free implementation of an FFT
  by Don Cross, available on the web at:

    http://www.intersrv.com/~dcross/fft.html

  The basic algorithm for his code was based on Numerical Recipes
  in Fortran.  I optimized his code further by reducing array
  accesses, caching the bit reversal table, and eliminating
  float-to-double conversions, and I added the routines to
  calculate a real FFT and a real power spectrum.

  Note: all of these routines use single-precision floats.
  I have found that in practice, floats work well until you
  get above 8192 samples.  If you need to do a larger FFT,
  you need to use doubles.

**********************************************************************/

#ifndef M_PI
#define	M_PI		3.14159265358979323846  /* pi */
#endif

#if (defined(__BCPLUSPLUS__) || defined(SOLARIS))
#define isgreater(x, y) (x > y)
#define isless(x, y)    (x < y)
#define fmax(x, y)      ((x >= y) ? x : y)
#define fmin(x, y)      ((x <= y) ? x : y)
#define fmaxf(x, y)     fmax(x, y)
#define fminf(x, y)     fmin(x, y)
#define fma(x, y, z)    (z + x * y)
#define fmaf(x, y, z)   fma(x, y, z)
#define fabsf(x)        fabs(x)
#define powf(x, y)      pow(x, y)
#define sqrtf(x)        sqrt(x)
#define expf(x)         exp(x)
#define logf(x)         log(x)
#define log10f(x)       log10(x)
#define sinf(x)         sin(x)
#define sinhf(x)        sinh(x)
#define cosf(x)         cos(x)
#define coshf(x)        cosh(x)
#define atanf(x)        atan(x)
#define atan2f(x, y)    atan2(x, y)
#define acosf(x)        acos(x)
#define acosh(x)        (log(x + sqrt(x * x - 1)))
#define acoshf(x)       (logf(x + sqrtf(x * x - 1)))
#define hypotf(x, y)    hypot(x, y)
#endif

#define DSP_MAXBESSEL       32L

static const int MaxFastBits = 16;

class TFFT {

static int **gFFTBitTable;
static int usetable;

public:

   TFFT(){
      InitFFT();
   };

   ~TFFT(){
      FreeFFT();
   };

   void InitFFT();
   int  FastReverseBits(int i, int NumBits);
     

/*
 * This is the function you will use the most often.
 * Given an array of floats, this will compute the power
 * spectrum by doing a Real FFT and then computing the
 * sum of the squares of the real and imaginary parts.
 * Note that the output array is half the length of the
 * input array, and that NumSamples must be a power of two.
 */

void PowerSpectrum(int NumSamples, double *In, double *Out);

/*
 * Computes an FFT when the input data is real but you still
 * want complex data as output.  The output arrays are half
 * the length of the input, and NumSamples must be a power of
 * two.
 */

void RealFFT(int NumSamples,
             double *RealIn, double *RealOut, double *ImagOut);

/*
 * Computes a FFT of complex input and returns complex output.
 * Currently this is the only function here that supports the
 * inverse transform as well.
 */

void FFT(int NumSamples,
         bool InverseTransform,
         double *RealIn, double *ImagIn, double *RealOut, double *ImagOut);

/*
 * Applies a windowing function to the data in place
 *
 * 0: Rectangular (no window)
 * 1: Bartlett    (triangular)
 * 2: Hamming
 * 3: Hanning
 * 4: Boghmann
 * 5: Pouasson
 */

void WindowFunc(int whichFunction, int NumSamples, double *data);

/*
 * Returns the name of the windowing function (for UI display)
 */

char *WindowFuncName(int whichFunction);

/*
 * Returns the number of windowing functions supported
 */

int NumWindowFuncs();
/*
* Free memory
*/
void FreeFFT();

};
#endif
