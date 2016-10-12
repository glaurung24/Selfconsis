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
#include "ChebyshevSolver.h"
#include "CPropertyExtractor.h"
#include "FileWriter.h"
#include "LDOS.h"
#include "Util.h"
#include <sstream>
#include "ParameterSet.h"
#include "FileParser.h"
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


    //Setup ChebyshevSolver
//    ChebyshevSolver cSolver;
//    cSolver.setModel(&model);
//    cSolver.setScaleFactor(SCALE_FACTOR);
//
//
//
//    //Create PropertyExtractor. The parameter are in order: The
//    //ChebyshevSolver, number of expansion coefficients used in the
//    //Cebyshev expansion, energy resolution with which the Green's function
//    // is evaluated, whether calculate expansion functions using a GPU or
//    //not, whether to evaluate the Green's function using a GPU or not,
//    //whether to use a lookup table for the Green's function or not
//    //(required if the Green's function is evaluated on a GPU), and the
//    //lower and upper bound between which the Green's function is evaluated
//    //(has to be inside the interval [-SCALE_FACTOR, SCALE_FACTOR]).
//    CPropertyExtractor pe(&cSolver,
//                NUM_COEFFICIENTS,
//                ENERGY_RESOLUTION,
//                false,
//                false,
//                true,
//                -1,
//                1);
//
//    //Extract local density of states and write to file
//    Property::LDOS *ldos = pe.calculateLDOS({IDX_X, SIZE_Y/2, IDX_SUM_ALL},
//                        {SIZE_X, 1, 2});
//
////	const int RANK = 1;
////	int dims[RANK] = {SIZE_X};
//
//
//    //Set filename and remove any file already in the folder
//    stringstream ss;
//    ss << "TBTKResults_" << counter_z << ".h5";
//    FileWriter::setFileName(ss.str());
//    FileWriter::clear();
//    FileWriter::writeLDOS(ldos);
//    delete ldos;
//        Util::Timer::tock();
	return 0;
}


