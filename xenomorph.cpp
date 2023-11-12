#include <iostream>
#include "xenomorph.h"
#include "tinyfft.h"

void XENOMORPH_EXPORT doMorfing( const double micData[ MIC_ARRAY_SIZE ][ BLOCK_SIZE ], double ret[ OUT_MAP_WIDTH ][ OUT_MAP_HEIGHT ] ) {

    /// TODO algo here
    for(int i=0; i<MIC_ARRAY_SIZE; i++)
        std::cout << i << "   " << COORD_X[i] << " <-> " << COORD_Y[i] << std::endl;
//    ret[42][42] = 4242.;
    for(int i=0; i<OUT_MAP_WIDTH; i++)
        for(int j=0; j<OUT_MAP_WIDTH; j++)
            ret[i][j] = 7*i+j;
}
