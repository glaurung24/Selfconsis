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
    Calculation::Init("input");
    Calculation::SetUpModel();
    Calculation::ScLoop(true);
    Calculation::CalcLDOS();

	return 0;
}


