#include <iostream>
#include <math.h>
#include "xenomorph.h"
#include "fft.h"

void XENOMORPH_EXPORT doMorfing( const double micData[ MIC_ARRAY_SIZE ][ BLOCK_SIZE ], double ret[ OUT_MAP_WIDTH ][ OUT_MAP_HEIGHT ] ) {

    double power_buff[BLOCK_SIZE];
    double real_fft[BLOCK_SIZE];
    double imag_fft[BLOCK_SIZE];
    double XH_real[MIC_ARRAY_SIZE];
    double XH_imag[MIC_ARRAY_SIZE];

    TFFT fft_obj;
    // find peaks
    fft_obj.PowerSpectrum(BLOCK_SIZE, (double *)micData[0], power_buff);
    int max_index = 0;
    double max_value = 0;
    for(int i=9; i<BLOCK_SIZE/2; i++){
        if(max_value<power_buff[i]){
            max_index = i;
            max_value = power_buff[i];
        }
    }
    double max_frequency = 1.0*max_index/BLOCK_SIZE*FD; // Hz
    std::cout << "F = " << max_frequency;

    // calc spectrums
    for(int i=0; i<MIC_ARRAY_SIZE; i++)
        fft_obj.RealFFT(BLOCK_SIZE, (double *)micData[i], real_fft, imag_fft);
    // find space peaks
    for(int i=0; i<OUT_MAP_WIDTH; i++)
        for(int j=0; j<OUT_MAP_WIDTH; j++){
            double phi = H_MIN + (H_MAX-H_MIN)*i/OUT_MAP_WIDTH;
            double theta = V_MIN + (V_MAX-V_MIN)*i/OUT_MAP_WIDTH;
            double real_sum = 0;
            double imag_sum = 0;
            for(int k=0; k<MIC_ARRAY_SIZE; k++){
                double delay = (COORD_X[k]*sin(phi) + COORD_Y[k]*sin(theta))/C;
                double phase = 2*M_PI*max_frequency*delay;
                real_sum += XH_real[k]*cos(phase) - XH_imag[k]*sin(phase);
                imag_sum += XH_real[k]*sin(phase) + XH_imag[k]*cos(phase);
                }
            ret[i][j] = sqrt( real_sum*real_sum + imag_sum*imag_sum );
            }
}
