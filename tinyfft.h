#include <complex>
#include "xenomorph_global.h"

//void XENOMORPH_EXPORT doMorfing( const double micData[ MIC_ARRAY_SIZE ][ BLOCK_SIZE ], double ret[ OUT_MAP_WIDTH ][ OUT_MAP_HEIGHT ] );

void tfft_init(const int k, double complex w[BLOCK_SIZE]);

void tfft_fft(int k, double complex A[], const double complex w[]);

void tfft_ifft(int k, double complex A[], const double complex w[]);

void tfft_convolver(int k, double complex A[], const double complex w[]);

//const int BLOCK_SIZE = 512;
