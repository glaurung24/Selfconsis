/** @package TBTKtemp
 *  @file main.cpp
 *  @brief Basic Chebyshev example
 *
 *  Basic example of using the Chebyshev method to solve a 2D tight-binding
 *  model with t = 1 and mu = -1. Lattice with edges and a size of 40x40 sites.
 *  Using 5000 Chebyshev coefficients and evaluating the Green's function with
 *  an energy resolution of 10000. Calculates LDOS at SIZE_X = 40 sites along
 *  the line y = SIZE_Y/2 = 20.
 *
 *  @author Kristofer Björnson
 */

#include <complex>
#include <iostream>
#include <memory>
#include "Model.h"
#include "ChebyshevSolver.h"
#include "CPropertyExtractor.h"
#include "FileWriter.h"
#include "LDOS.h"

using namespace std;
using namespace TBTK;

const complex<double> i(0, 1);

int main(int argc, char **argv){
	//Lattice size
	const int N = 10;
	const int SIZE_X = 4*N;
	const int SIZE_Y = 2*N+1;
	const int SPIN_D = 4;

    //Spin index:
    // 0: particle up
    // 1: particle down
    // 2: Hole up
    // 3: Hole down


	//Model parameters.
	complex<double> mu = -1.0;
	complex<double> t = 1.0;
	complex<double> z = 2.0; //Zeeman coupling
	complex<double> delta[SIZE_X][SIZE_Y];
	complex<double> deltaStart = 0.3;

	bool isMagnetized[SIZE_X][SIZE_Y];



	for(int i = 0; i < SIZE_X; i++){
        for(int j = 0; j < SIZE_Y; j++){
            delta[i][j]=deltaStart;
            isMagnetized[i][j]= false;
        }
	}

	for(int x = SIZE_X/4; x < 3*SIZE_X/4; x++){
        isMagnetized[x][SIZE_Y/2]=true;
	}

	//Create model and set up hopping parameters
	Model model;
	for(int x = 0; x < SIZE_X; x++){
		for(int y = 0; y < SIZE_Y; y++){
			for(int s = 0; s < SPIN_D/2; s++){

//------------------------chemical Potential-----------------------------------
				//Add hopping amplitudes corresponding to chemical potential
				model.addHA(HoppingAmplitude(-mu,	{x, y, s},	{x, y, s}));
				model.addHA(HoppingAmplitude(mu,	{x, y, s+2},	{x, y, s+2}));

//-------------------interaction term------------------------------------------

                model.addHAAndHC(HoppingAmplitude(delta[x][y]*2.0*(0.5-s), {x,y,s}, {x,y,(3-s)}));


//------------------------Rashba coupling or hopping term--------------------------------------
				//Add hopping parameters corresponding to t
				if(x+1 < SIZE_X){
					model.addHAAndHC(HoppingAmplitude(-t,	{(x+1)%SIZE_X, y, s},	{x, y, s}));
					model.addHAAndHC(HoppingAmplitude(t,	{x, y, s+2},{(x+1)%SIZE_X, y, s+2}));
//					model.addHAAndHC(HoppingAmplitude(t, {(x+1)%SIZE_X, y, s+2}, {x, y, s+2})); //same results with this line as with the line above
                }
				if(y+1 < SIZE_Y){
					model.addHAAndHC(HoppingAmplitude(-t,	{x, (y+1)%SIZE_Y, s},	{x, y, s}));
					model.addHAAndHC(HoppingAmplitude(t,  {x, y, s+2}, {x, (y+1)%SIZE_Y, s+2}));
//					model.addHAAndHC(HoppingAmplitude(t,   {x, (y+1)%SIZE_Y, s+2},{x, y, s+2}));
                }

//---------------------------Zeeman term------------------------------------------
                if(isMagnetized[x][y]){
                    model.addHA(HoppingAmplitude(z*2.0*(0.5-s), {x, y, s}, {x, y, s}));
                    model.addHA(HoppingAmplitude(-z*2.0*(0.5-s), {x, y, s+2}, {x, y, s+2}));
                }


			}
		}
	}

	//Construct model
	model.construct();

	//TODO ask about fermi statistic

	//Chebyshev expansion parameters.
	const int NUM_COEFFICIENTS = 1000;
	const int ENERGY_RESOLUTION = 2000;
	const double SCALE_FACTOR = 10.;

	//Setup ChebyshevSolver
	ChebyshevSolver cSolver;
	cSolver.setModel(&model);
	cSolver.setScaleFactor(SCALE_FACTOR);

	//Set filename and remove any file already in the folder
	FileWriter::setFileName("TBTKResults.h5");
	FileWriter::clear();

	//Create PropertyExtractor. The parameter are in order: The
	//ChebyshevSolver, number of expansion coefficients used in the
	//Cebyshev expansion, energy resolution with which the Green's function
	// is evaluated, whether calculate expansion functions using a GPU or
	//not, whether to evaluate the Green's function using a GPU or not,
	//whether to use a lookup table for the Green's function or not
	//(required if the Green's function is evaluated on a GPU), and the
	//lower and upper bound between which the Green's function is evaluated
	//(has to be inside the interval [-SCALE_FACTOR, SCALE_FACTOR]).
	CPropertyExtractor pe(&cSolver,
				NUM_COEFFICIENTS,
				ENERGY_RESOLUTION,
				false,
				false,
				true,
				-1,
				1);

	//Extract local density of states and write to file
	Property::LDOS *ldos = pe.calculateLDOS({IDX_X, SIZE_Y/2, IDX_SUM_ALL},
						{SIZE_X, 1, 2});
//	const int RANK = 1;
//	int dims[RANK] = {SIZE_X};
	FileWriter::writeLDOS(ldos);
	delete ldos;

	return 0;
}
