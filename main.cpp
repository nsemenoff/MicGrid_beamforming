#include <iostream>
#include <math.h>
#include "xenomorph.h"
#include "fft.h"

#define FREQ 1000

int main() {

    double mphData[ MIC_ARRAY_SIZE ][ BLOCK_SIZE ]{{0.,},}; //
    double result[ OUT_MAP_WIDTH ][ OUT_MAP_HEIGHT ]{{0.,},}; // начало координат [0][0] -- верхний левый угол

    for(int i=0; i<MIC_ARRAY_SIZE; i++)
        for(int j=0; j<BLOCK_SIZE; j++)
            mphData[i][j] = 25000*sin(2*M_PI*FREQ*j/FD);

    doMorfing( mphData, result );

    std::cout << result[12][12];

    return 0;
}
