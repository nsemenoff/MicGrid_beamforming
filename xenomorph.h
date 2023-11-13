#ifndef XENOMORPH_H
#define XENOMORPH_H

#include "xenomorph_global.h"

#ifndef MIC_ARRAY_SIZE
#define MIC_ARRAY_SIZE 49
#endif

const int BLOCK_SIZE = 512;
const int OUT_MAP_WIDTH = 16;
const int OUT_MAP_HEIGHT = 16;
const double H_MIN = -60;
const double H_MAX = 60;
const double V_MIN = -60;
const double V_MAX = 60;

const double FD = 48000;
const double C = 335;
const double COORD_X[MIC_ARRAY_SIZE] = {-0.22, -0.12, -0.05, 0.00, 0.05, 0.12, 0.22,
                            -0.22, -0.12, -0.05, 0.00, 0.05, 0.12, 0.22,
                            -0.22, -0.12, -0.05, 0.00, 0.05, 0.12, 0.22,
                            -0.22, -0.12, -0.05, 0.00, 0.05, 0.12, 0.22,
                            -0.22, -0.12, -0.05, 0.00, 0.05, 0.12, 0.22,
                            -0.22, -0.12, -0.05, 0.00, 0.05, 0.12, 0.22,
                            -0.22, -0.12, -0.05, 0.00, 0.05, 0.12, 0.22};
const double COORD_Y[MIC_ARRAY_SIZE] = {-0.22, -0.22, -0.22, -0.22, -0.22, -0.22, -0.22,
                            -0.12, -0.12, -0.12, -0.12, -0.12, -0.12, -0.12,
                            -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05,
                             0.00,  0.00,  0.00,  0.00,  0.00,  0.00,  0.00,
                             0.05,  0.05,  0.05,  0.05,  0.05,  0.05,  0.05,
                             0.12,  0.12,  0.12,  0.12,  0.12,  0.12,  0.12,
                             0.22,  0.22,  0.22,  0.22,  0.22,  0.22,  0.22};

void XENOMORPH_EXPORT doMorfing( const double micData[ MIC_ARRAY_SIZE ][ BLOCK_SIZE ], double ret[ OUT_MAP_WIDTH ][ OUT_MAP_HEIGHT ] );

#endif // XENOMORPH_H
