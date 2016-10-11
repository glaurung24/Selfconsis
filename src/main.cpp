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

using namespace std;
using namespace TBTK;

const complex<double> i(0, 1);





//complex<double>*** delta;
//complex<double> ***delta_current, ***delta_old, ***delta_tmp;


complex<double> funcDelta(Index, Index);
void initDelta();
void scLoop();

int main(int argc, char **argv){
	//Lattice size
	const int N = 20;
    int SIZE_X = 4*N;
    int SIZE_Y = 2*N+1;
    int SPIN_D = 4;

    //Spin index:
    // 0: particle up
    // 1: particle down
    // 2: Hole up
    // 3: Hole down

    Util::ParameterSet* ps = FileParser::readParameterSet("input");
    int counter_z = ps->getInt("counter_z");

    //Model parameters.
//        Util::Timer::tick("setUpModel");

    complex<double> mu = -4.0;
    complex<double> t = 1.0;
    complex<double> z = ps->getComplex("z"); //Zeeman coupling
    complex<double> delta[SIZE_X][SIZE_Y];
    complex<double> deltaStart = 0.3;
    complex<double> alpha=0.3;

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

//-------------------BCS interaction term------------------------------------------

//                model.addHAAndHC(HoppingAmplitude(delta[x][y]*2.0*(0.5-s), {x,y,s}, {x,y,(3-s)}));
                model.addHAAndHC(HoppingAmplitude(funcDelta, {x,y,s}, {x,y,(3-s)}));


//------------------------Nearest neighbour hopping term--------------------------------------
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

//------------------------Rashba hopping term--------------------------------------

                if(x+1 < SIZE_X){
    //                    model.addHAAndHC(HoppingAmplitude(alpha*2.0*(0.5-s), {(x+1)%SIZE_X, y, s*2},	{x, y, s*2+1}));
    //                    model.addHAAndHC(HoppingAmplitude(-alpha*2.0*(0.5-s), {x, y, s*2},	{(x+1)%SIZE_X, y, s*2+1}));
                    model.addHAAndHC(HoppingAmplitude(alpha *2.0*(0.5-s), {(x+1)%SIZE_X,y,(s+1)%2}, {x,y,s}));
                    model.addHAAndHC(HoppingAmplitude(-alpha *2.0*(0.5-s), {x,y,s+2}, {(x+1)%SIZE_X,y,(s+1)%2+2}));
                }

                if(y+1 < SIZE_Y){
    //                    model.addHAAndHC(HoppingAmplitude(i*alpha*2.0*(0.5-s),	{x, (y+1)%SIZE_Y, s*2},	{x, y, s*2+1}));
    //                    model.addHAAndHC(HoppingAmplitude(-i*alpha*2.0*(0.5-s),  {x, y, s*2}, {x, (y+1)%SIZE_Y, s*2+1}));
                    model.addHAAndHC(HoppingAmplitude(i*alpha, {x,(y+1)%SIZE_Y,(s+1)%2}, {x,y,s}));
                    model.addHAAndHC(HoppingAmplitude(-i*alpha, {x,y,s+2}, {x,(y+1)%SIZE_Y,(s+1)%2+2}));

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
//        Util::Timer::tock();

    //TODO ask about fermi statistic
//        Util::Timer::tick("Calculation");
    //Chebyshev expansion parameters.
    const int NUM_COEFFICIENTS = 1000;
    const int ENERGY_RESOLUTION = 2000;
    const double SCALE_FACTOR = 10.;

    //Setup ChebyshevSolver
    ChebyshevSolver cSolver;
    cSolver.setModel(&model);
    cSolver.setScaleFactor(SCALE_FACTOR);



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


    //Set filename and remove any file already in the folder
    stringstream ss;
    ss << "TBTKResults_" << counter_z << ".h5";
    FileWriter::setFileName(ss.str());
    FileWriter::clear();
    FileWriter::writeLDOS(ldos);
    delete ldos;
//        Util::Timer::tock();
	return 0;
}


complex<double>funcDelta(Index to, Index from)
{
    int from_x = from.at(0);
    int from_y = from.at(1);
    int from_s = from.at(2);
//    int to_x = to.at(0);
//    int to_y = to.at(1);
//    int to_s = to.at(2);


}

void scLoop()
{
//    #pragma omp parallel for
//    complex<double> pairFunction = pe.calculateExpectationValue({x,y,3},{x,y,0});

//    delta_new[x][y] = -interactionValueBCS*pairFunction;
}

void initDelta()
{


}
