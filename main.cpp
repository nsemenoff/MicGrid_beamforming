#include <iostream>
#include "xenomorph.h"

int main() {

    double mphData[ MIC_ARRAY_SIZE ][ BLOCK_SIZE ]{{0.,},}; // 
    double result[ OUT_MAP_WIDTH ][ OUT_MAP_HEIGHT ]{{0.,},}; // начало координат [0][0] -- верхний левый угол

    doMorfing( mphData, result );

    
    std::cout << result[12][12];

    return 0;
}
