/**********************************************************************

  FFT.cpp

  Dominic Mazzoni

  September 2000

  This file contains a few FFT routines, including a real-FFT
  routine that is almost twice as fast as a normal complex FFT,
  and a power spectrum routine when you know you don't care
  about phase information.

  Some of this code was based on a free implementation of an FFT
  by Don Cross, available on the web at:

    http://www.intersrv.com/~dcross/fft.html

  The basic algorithm for his code was based on Numerican Recipes
  in Fortran.  I optimized his code further by reducing array
  accesses, caching the bit reversal table, and eliminating
  float-to-double conversions, and I added the routines to
  calculate a real FFT and a real power spectrum.

**********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <vector>

using namespace std;

#include "fft.h"

int **TFFT::gFFTBitTable = NULL;
int TFFT::usetable = 0;

static int IsPowerOfTwo(int x)
{
   if (x < 2)
      return false;

   if (x & (x - 1))             /* Thanks to 'byang' for this cute trick! */
      return false;

   return true;
}

static int NumberOfBitsNeeded(int PowerOfTwo)
{
   int i;

   if (PowerOfTwo < 2) {
      fprintf(stderr, "Error: FFT called with size %d\n", PowerOfTwo);
      exit(1);
   }

   for (i = 0;; i++)
      if (PowerOfTwo & (1 << i))
         return i;
}

static int ReverseBits(int index, int NumBits)
{
   int i, rev;

   for (i = rev = 0; i < NumBits; i++) {
      rev = (rev << 1) | (index & 1);
      index >>= 1;
   }

   return rev;
}

void TFFT::InitFFT()
{
   if(usetable == 0){
      gFFTBitTable = new int *[MaxFastBits];
      int len = 2;
      for (int b = 1; b <= MaxFastBits; b++) {
         gFFTBitTable[b - 1] = new int[len];
         for (int i = 0; i < len; i++)
            gFFTBitTable[b - 1][i] = ReverseBits(i, b);
         len <<= 1;
      }
   }// if
   usetable++;
}

void TFFT::FreeFFT()
{
   usetable--;
   if(!usetable && gFFTBitTable){
      for (int b = 1; b <= MaxFastBits; b++)
         delete[] gFFTBitTable[b - 1];
      delete[] gFFTBitTable;
      gFFTBitTable = NULL;
   }// if
}

int TFFT::FastReverseBits(int i, int NumBits)
{
   if (NumBits <= MaxFastBits)
      return gFFTBitTable[NumBits - 1][i];
   else
      return ReverseBits(i, NumBits);
}

/*
 * Complex Fast Fourier Transform
 */

void TFFT::FFT(int NumSamples,
         bool InverseTransform,
         double *RealIn, double *ImagIn, double *RealOut, double *ImagOut)
{
   int NumBits;                 /* Number of bits needed to store indices */
   int i, j, k, n;
   int BlockSize, BlockEnd;

   double angle_numerator = 2.0 * M_PI;
   double tr, ti;                /* temp real, temp imaginary */

   if (!IsPowerOfTwo(NumSamples)) {
      fprintf(stderr, "%d is not a power of two\n", NumSamples);
      exit(1);
   }

   if (!gFFTBitTable)
      InitFFT();

   if (InverseTransform)
      angle_numerator = -angle_numerator;

   NumBits = NumberOfBitsNeeded(NumSamples);

   /*
    **   Do simultaneous data copy and bit-reversal ordering into outputs...
    */

   for (i = 0; i < NumSamples; i++) {
      j = FastReverseBits(i, NumBits);
      RealOut[j] = RealIn[i];
      ImagOut[j] = (ImagIn == NULL) ? 0.0 : ImagIn[i];
   }

   /*
    **   Do the FFT itself...
    */

   BlockEnd = 1;
   for (BlockSize = 2; BlockSize <= NumSamples; BlockSize <<= 1) {

      double delta_angle = angle_numerator / (double) BlockSize;

      double sm2 = sin(-2 * delta_angle);
      double sm1 = sin(-delta_angle);
      double cm2 = cos(-2 * delta_angle);
      double cm1 = cos(-delta_angle);
      double w = 2 * cm1;
      double ar0, ar1, ar2, ai0, ai1, ai2;

      for (i = 0; i < NumSamples; i += BlockSize) {
         ar2 = cm2;
         ar1 = cm1;

         ai2 = sm2;
         ai1 = sm1;

         for (j = i, n = 0; n < BlockEnd; j++, n++) {
            ar0 = w * ar1 - ar2;
            ar2 = ar1;
            ar1 = ar0;

            ai0 = w * ai1 - ai2;
            ai2 = ai1;
            ai1 = ai0;

            k = j + BlockEnd;
            tr = ar0 * RealOut[k] - ai0 * ImagOut[k];
            ti = ar0 * ImagOut[k] + ai0 * RealOut[k];

            RealOut[k] = RealOut[j] - tr;
            ImagOut[k] = ImagOut[j] - ti;

            RealOut[j] += tr;
            ImagOut[j] += ti;
         }
      }

      BlockEnd = BlockSize;
   }

   /*
      **   Need to normalize if inverse transform...
    */

   if (InverseTransform) {
      double denom = (double) NumSamples;

      for (i = 0; i < NumSamples; i++) {
         RealOut[i] /= denom;
         ImagOut[i] /= denom;
      }
   }
}

/*
 * Real Fast Fourier Transform
 *
 * This function was based on the code in Numerical Recipes in C.
 * In Num. Rec., the inner loop is based on a single 1-based array
 * of interleaved real and imaginary numbers.  Because we have two
 * separate zero-based arrays, our indices are quite different.
 * Here is the correspondence between Num. Rec. indices and our indices:
 *
 * i1  <->  real[i]
 * i2  <->  imag[i]
 * i3  <->  real[n/2-i]
 * i4  <->  imag[n/2-i]
 */

void TFFT::RealFFT(int NumSamples, double *RealIn, double *RealOut, double *ImagOut)
{
   int Half = NumSamples / 2;
   int i;

   double theta = M_PI / Half;

   vector<double> tmpReal(Half);
   vector<double> tmpImag(Half);

   
   for (i = 0; i < Half; i++) {
      tmpReal[i] = RealIn[2 * i];
      tmpImag[i] = RealIn[2 * i + 1];
   }

   FFT(Half, 0, &tmpReal[0], &tmpImag[0], RealOut, ImagOut);

   double wtemp = double (sin(0.5 * theta));

   double wpr = -2.0 * wtemp * wtemp;
   double wpi = double (sin(theta));
   double wr = 1.0 + wpr;
   double wi = wpi;

   int i3;

   double h1r, h1i, h2r, h2i;

   for (i = 1; i < Half / 2; i++) {

      i3 = Half - i;

      h1r = 0.5 * (RealOut[i] + RealOut[i3]);
      h1i = 0.5 * (ImagOut[i] - ImagOut[i3]);
      h2r = 0.5 * (ImagOut[i] + ImagOut[i3]);
      h2i = -0.5 * (RealOut[i] - RealOut[i3]);

      RealOut[i] = h1r + wr * h2r - wi * h2i;
      ImagOut[i] = h1i + wr * h2i + wi * h2r;
      RealOut[i3] = h1r - wr * h2r + wi * h2i;
      ImagOut[i3] = -h1i + wr * h2i + wi * h2r;

      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
   }

   RealOut[0] = (h1r = RealOut[0]) + ImagOut[0];
   ImagOut[0] = h1r - ImagOut[0];

}

/*
 * PowerSpectrum
 *
 * This function computes the same as RealFFT, above, but
 * adds the squares of the real and imaginary part of each
 * coefficient, extracting the power and throwing away the
 * phase.
 *
 * For speed, it does not call RealFFT, but duplicates some
 * of its code.
 */

void TFFT::PowerSpectrum(int NumSamples, double *In, double *Out)
{
   int Half = NumSamples / 2;
   int i;

   double theta = M_PI / Half;

   vector<double> tmpReal(Half);
   vector<double> tmpImag(Half);
   vector<double> RealOut(Half);
   vector<double> ImagOut(Half);
   
   for (i = 0; i < Half; i++) {
      tmpReal.at(i) = In[2 * i];
      tmpImag.at(i) = In[2 * i + 1];
   }

   FFT(Half, 0, &tmpReal[0], &tmpImag[0], &RealOut[0], &ImagOut[0]);

   double wtemp = double (sin(0.5 * theta));

   double wpr = -2.0 * wtemp * wtemp;
   double wpi = double (sin(theta));
   double wr = 1.0 + wpr;
   double wi = wpi;

   int i3;

   double h1r, h1i, h2r, h2i, rt, it;

   for (i = 1; i < Half / 2; i++) {

      i3 = Half - i;

      h1r = 0.5 * (RealOut.at(i) + RealOut[i3]);
      h1i = 0.5 * (ImagOut.at(i) - ImagOut[i3]);
      h2r = 0.5 * (ImagOut.at(i) + ImagOut[i3]);
      h2i = -0.5 * (RealOut.at(i) - RealOut[i3]);

      rt = h1r + wr * h2r - wi * h2i;
      it = h1i + wr * h2i + wi * h2r;

      Out[i] = rt * rt + it * it;

      rt = h1r - wr * h2r + wi * h2i;
      it = -h1i + wr * h2i + wi * h2r;

      Out[i3] = rt * rt + it * it;

      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
   }

   rt = (h1r = RealOut[0]) + ImagOut[0];
   it = h1r - ImagOut[0];
   Out[0] = rt * rt + it * it;

   rt = RealOut[Half / 2];
   it = ImagOut[Half / 2];
   Out[Half / 2] = rt * rt + it * it;

//   delete[]tmpReal;
//   delete[]tmpImag;
//   delete[]RealOut;
//   delete[]ImagOut;
}

/*
 * Windowing Functions
 */

int TFFT::NumWindowFuncs()
{
   return 6;
}

double ChebyshevPolynom (double dOrder, double dValue)
{
    if (fabs(dValue) <= 1.0)
    {
        return cos(dOrder * acos(dValue));
    }
    else
    {
        return cosh(dOrder * acos(dValue));
    }
}

void MinMax (double *dpMin, double *dpMax, const double *dpSrc,
    long lCount)
{
    long lLoopCntr;
    double dTempVal;
    double dTempMin = 1;
    double dTempMax = 0;

    for (lLoopCntr = 0L; lLoopCntr < lCount; lLoopCntr++)
    {
        dTempVal = dpSrc[lLoopCntr];
        #ifndef _ISOC9X_SOURCE
        if (dTempVal < dTempMin)
        {
            dTempMin = dTempVal;
        }
        if (dTempVal > dTempMax)
        {
            dTempMax = dTempVal;
        }
        #else
        dTempMin = fmin(dTempVal, dTempMin);
        dTempMax = fmax(dTempVal, dTempMax);
        #endif
    }
    *dpMin = dTempMin;
    *dpMax = dTempMax;
}

double Multiple (long lValue)
{
    long lLoopCntr;
    double dMult = 1.0;

    for (lLoopCntr = 1L; lLoopCntr <= lValue; lLoopCntr++)
    {
        dMult *= (double) lLoopCntr;
    }
    return dMult;
}

double ModZeroBessel (double dValue)
{
    long lLoopCntr;
    double dMZBessel = 0.0;
    double dHalfValue;

    dHalfValue = dValue / 2.0;
    for (lLoopCntr = 0L; lLoopCntr <= DSP_MAXBESSEL; lLoopCntr++)
    {
        dMZBessel += pow(
            pow(dHalfValue, lLoopCntr) / Multiple(lLoopCntr), 2.0);
    }
    return dMZBessel;
}

char *TFFT::WindowFuncName(int whichFunction)
{
   switch (whichFunction) {
   default:
   case 0:
      return "Rectangular";
   case 1:
      return "Bartlett";
   case 2:
      return "Hamming";
   case 3:
      return "Hanning";
   case 4:
      return "Boghmann";
   case 5:
      return "Puasson";
   }
}

void TFFT::WindowFunc(int whichFunction, int NumSamples, double *in)
{
   int i;
   double N  = (double)NumSamples, 
          N2 = (double)(NumSamples / 2);

   if (whichFunction == 5) {
      // Pouasson
      double K = 0.1;
      for (i = 0; i < NumSamples / 2; i++) {
         in[i                   ] *= 1 - pow(K , 2 * i / N);
         in[i + (NumSamples / 2)] *= 1 - pow(K , 2 * (N - i) / N);
      }
   }

   if (whichFunction == 4) {
      // Boghmann
/*
      double w, ws, OneDivPI = (double)1.0 / M_PI;
      for (i = 0; i < NumSamples / 2; i++) {
         w = 2 * M_PI * i / N;
         in[i                   ] *= (2 * i / N) * cos(w) + OneDivPI * sin(w * i);
         ws = w * (i - N2);
         in[i + (NumSamples / 2)] *= (2 - 2 * i / N) * cos(ws) + OneDivPI * sin(w);
      }
*/

{
    long lLoopCntr, lSize = NumSamples;

    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        in[lLoopCntr] *= 0.42323 -
            0.49855 *
            cos((2.0 * M_PI * (double) lLoopCntr) / (double) lSize) +
            0.07922 *
            cos((4.0 * M_PI * (double) lLoopCntr) / (double) lSize);
    }
}

/*
{
    long lLoopCntr, lSize = NumSamples;
    double dA0;
    double dA1;
    double dA2;

    dA0 = 7938.0 / 18608.0;
    dA1 = 9240.0 / 18608.0;
    dA2 = 1430.0 / 18608.0;
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        in[lLoopCntr] *= dA0 -
            dA1 * cos((2.0 * M_PI * (double) lLoopCntr) / (double) lSize) +
            dA2 * cos((4.0 * M_PI * (double) lLoopCntr) / (double) lSize);
    }
}
*/
/*
{
    long lLoopCntr, lSize = NumSamples;
    float fHalfN;
    float fPiAlpha, fAlpha = 1.0;

    fHalfN = (float) lSize / 2.0F;
    fPiAlpha = M_PI * fAlpha;
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        in[lLoopCntr] *=
            ModZeroBessel(fPiAlpha * 
            sqrt(1.0F - pow(((float) lLoopCntr - fHalfN) / fHalfN, 2.0F))) /
            ModZeroBessel(fPiAlpha);
    }
}
*/
/*
{
    long lLoopCntr, lSize = NumSamples;
    long lHalfSize;

    lHalfSize = lSize / 2L;
    for (lLoopCntr = 0L; lLoopCntr < lSize; lLoopCntr++)
    {
        in[lLoopCntr] *= 0.5 *
            (1.0 + cos((float) (lLoopCntr - lHalfSize) * M_PI / lHalfSize));
    }
}
*/

   }// if

   if (whichFunction == 1) {
      // Bartlett (triangular) window
      for (i = 0; i < NumSamples / 2; i++) {
         double v = 2.0 * i / N2;
         in[i                   ] *= v;
         in[i + (NumSamples / 2)] *= 2.0 - v;
      }
   }

   if (whichFunction == 2) {
      // Hamming
      for (i = 0; i < NumSamples; i++)
         in[i] *= 0.54 - 0.46 * cos(2 * M_PI * i / (NumSamples - 1));
   }

   if (whichFunction == 3) {
      // Hanning
      for (i = 0; i < NumSamples; i++)
         in[i] *= 0.50 - 0.50 * cos(2 * M_PI * i / (NumSamples - 1));
   }
}
