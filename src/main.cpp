/** @package TBTKtemp
 *  @file main.cpp
 *  @brief Basic Chebyshev example
 *
 *  TODO
 *
 *
 *  @author template by  Kristofer Björnson
 */

#include <complex>
#include <iostream>
#include <memory>
#include "Model.h"
#include "Util.h"
#include <sstream>

#include "Calculation.h"

using namespace std;
using namespace TBTK;







//complex<double>*** delta;
//complex<double> ***delta_current, ***delta_old, ***delta_tmp;


complex<double> funcDelta(Index, Index);
void initDelta();
void scLoop();

int main(int argc, char **argv){
	//Lattice size
    Calculation::Init();
    Calculation::SetUpModel();
    Calculation::ScLoop();
//    Calculation::Delete(); //TODO uncomment

    //Spin index:
    // 0: particle up
    // 1: particle down
    // 2: Hole up
    // 3: Hole down



    //Model parameters.
//        Util::Timer::tick("setUpModel");




//        Util::Timer::tock();

    //TODO ask about fermi statistic
//        Util::Timer::tick("Calculation");
    //Chebyshev expansion parameters.



//        Util::Timer::tock();
	return 0;
}


