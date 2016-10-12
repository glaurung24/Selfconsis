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


int Calculation::N;
int Calculation::SIZE_X;
int Calculation::SIZE_Y;
int Calculation::SPIN_D;

complex<double> Calculation::mu;
complex<double> Calculation::t;
complex<double> Calculation::z; //Zeeman coupling
complex<double> Calculation::deltaStart;
complex<double> Calculation::alpha;
complex<double> Calculation::couplingPotential;

Model Calculation::model;

int Calculation::NUM_COEFFICIENTS;
int Calculation::ENERGY_RESOLUTION;
double Calculation::SCALE_FACTOR;
complex<double>** Calculation::deltaCurrent;
complex<double>** Calculation::deltaOld;
bool Calculation::checkInit = false;
int Calculation::numberSCRuns;
int Calculation::epsDelta;
bool Calculation::verbose;


void Calculation::Init()
{
        checkInit = true;
        N = 2;
        SIZE_X = 4*N;
        SIZE_Y = 2*N+1;
        SPIN_D = 4;

        mu = -4.0;
        t = 1.0;
        z = 0.5; //Zeeman coupling

        deltaStart = 0.3;
        alpha = 0.3;
        couplingPotential = 0.1;

        InitDelta();

        epsDelta = 0.001;
        numberSCRuns = 1;
        verbose = false;
}

void Calculation::InitDelta()
{
    deltaCurrent = new complex<double>*[SIZE_Y];
    deltaOld = new complex<double>*[SIZE_Y];

    int j=0;
    for(int i=0; i < SIZE_X; i++)
    {
        deltaCurrent[i] = new complex<double>[SIZE_X];
        deltaOld[i] = new complex<double>[SIZE_X];
        deltaCurrent[i][j] = deltaStart;
        deltaOld[i][j] = deltaStart;
        j++;
    }
}

void Calculation::Delete() //TODO seg fault!!!
{
    for(int i=0; i < SIZE_Y; i++)
    {
        delete deltaCurrent[i];
        delete deltaOld[i];
    }
    delete deltaCurrent;
    delete deltaOld;
}


void Calculation::Init(std::string input_file) //TODO
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
    //dtor //TODO destroy delta
}


void Calculation::SetUpModel()
{
    if(checkInit)
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
                    model.addHAAndHC(HoppingAmplitude(&Calculation::FuncDelta, {x,y,s}, {x,y,(3-s)}));


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
    }
    else
    {
        throw runtime_error("Calculation was not Initialised. Run Calculation::Init() before setting up Model");
    }
    //Construct model
    model.construct();
}


complex<double> Calculation::FuncDelta(Index to, Index from)
{
    int from_x = from.at(0);
    int from_y = from.at(1);
    int from_s = from.at(2);
    return deltaCurrent[from_x][from_y]*2.0*(0.5-from_s);
//    int to_x = to.at(0);
//    int to_y = to.at(1);
//    int to_s = to.at(2);


}



complex<double> DummyFunc(Index, Index)
{
    return 0.1+0.1*i;
}


void Calculation::ScLoop()
{
//    #pragma omp parallel for
//    complex<double> pairFunction = pe.calculateExpectationValue({x,y,3},{x,y,0});

//    delta_new[x][y] = -interactionValueBCS*pairFunction;

    ChebyshevSolver cSolver;
    cSolver.setModel(&model);
    cSolver.setScaleFactor(SCALE_FACTOR);

    cout << deltaCurrent[5][0] << endl;


    CPropertyExtractor pe(&cSolver,
                NUM_COEFFICIENTS,
                ENERGY_RESOLUTION,
                false,
                false,
                true,
                -SCALE_FACTOR,
                SCALE_FACTOR);

    int loopCounter = 0;
    while(loopCounter < numberSCRuns)
    {
        //TODO swap delta
        for(int x=0; x < SIZE_X; x++)
        {
            for(int y=0; y < SIZE_Y; y++)
            {
//                pe.calculateExpectationValue({x,y,3},{x,y,0});
                cout << x << ", " << y << endl;
                deltaCurrent[x][y] = DummyFunc({x,y,3},{x,y,0})*couplingPotential;
                cout << deltaCurrent[x][y] << endl;
            }
        }
        cout << "bla" << endl;
        if(!verbose)
        {
            cout << "Loop number: " << loopCounter << endl;
        }
        loopCounter++;
    }
}


