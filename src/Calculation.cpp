#include "Calculation.h"
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


Calculation::Calculation()
{
        N = 20;
        SIZE_X = 4*N;
        SIZE_Y = 2*N+1;
        SPIN_D = 4;

        mu = -4.0;
        t = 1.0;
        z = 0.5; //Zeeman coupling

        deltaStart = 0.3;
        alpha=0.3;

        InitDelta();
}

void Calculation::InitDelta()
{
        delta = new complex<double>**[2];
        delta[0] = new complex<double>*[SIZE_X];
        delta[1] = new complex<double>*[SIZE_X];

        int j=0;
        for(int i=0; i < SIZE_Y; i++)
        {
            delta[0][i] = new complex<double>[SIZE_Y];
            delta[1][i] = new complex<double>[SIZE_Y];
            delta[0][i][j] = deltaStart;
            j++;
        }
        *delta_current = delta[0];
}


Calculation::Calculation(std::string input_file)
{

    Util::ParameterSet* ps = FileParser::readParameterSet("input");
//    counter_z = ps->getInt("counter_z");
    N = 20;
    SIZE_X = 4*N;
    SIZE_Y = 2*N+1;
    SPIN_D = 4;

    mu = -4.0;
    t = 1.0;
    z = ps->getComplex("z"); //Zeeman coupling
    deltaStart = 0.3;
    alpha=0.3;

    InitDelta();
}


Calculation::~Calculation()
{
    //dtor
}


void Calculation::SetUpModel()
{
    bool isMagnetized[SIZE_X][SIZE_Y];



    for(int x = SIZE_X/4; x < 3*SIZE_X/4; x++){
        isMagnetized[x][SIZE_Y/2]=true;
    }

    //Create model and set up hopping parameters

    for(int x = 0; x < SIZE_X; x++){
        for(int y = 0; y < SIZE_Y; y++){
            for(int s = 0; s < SPIN_D/2; s++){

//------------------------chemical Potential-----------------------------------
                //Add hopping amplitudes corresponding to chemical potential
                model.addHA(HoppingAmplitude(-mu,	{x, y, s},	{x, y, s}));
                model.addHA(HoppingAmplitude(mu,	{x, y, s+2},	{x, y, s+2}));

//-------------------BCS interaction term------------------------------------------

//                model.addHAAndHC(HoppingAmplitude(delta[x][y]*2.0*(0.5-s), {x,y,s}, {x,y,(3-s)}));
                model.addHAAndHC(HoppingAmplitude(this->funcDelta, {x,y,s}, {x,y,(3-s)}));


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
}


complex<double> Calculation::funcDelta(Index to, Index from)
{
    int from_x = from.at(0);
    int from_y = from.at(1);
    int from_s = from.at(2);
    return *delta_current[from_x][from_y]*2.0*(0.5-from_s);
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

