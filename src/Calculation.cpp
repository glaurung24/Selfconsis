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
complex<double>** Calculation::deltaNew;
complex<double>** Calculation::deltaOld;
bool Calculation::checkInit = false;
bool Calculation::modelSetUp = false;
int Calculation::numberSCRuns;
int Calculation::epsDelta;
bool Calculation::verbose;

unique_ptr<ChebyshevSolver> Calculation::cSolver;
unique_ptr<CPropertyExtractor> Calculation::pe;

string Calculation::fileName;


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
        couplingPotential = 1;

        InitDelta();

        epsDelta = 0.001;
        numberSCRuns = 5;
        verbose = false;

        cSolver = nullptr;
        pe = nullptr;

        NUM_COEFFICIENTS = 1000;
        ENERGY_RESOLUTION = 2000;
        SCALE_FACTOR = 10.;

        fileName = "TBTKResults.h5";
}

void Calculation::InitDelta()
{
    deltaNew = new complex<double>*[SIZE_X];
    deltaOld = new complex<double>*[SIZE_X];

    for(int i=0; i < SIZE_X; i++)
    {
        deltaNew[i] = new complex<double>[SIZE_Y];
        deltaOld[i] = new complex<double>[SIZE_Y];
        for(int j=0; j < SIZE_Y; j++)
        {
            deltaNew[i][j] = deltaStart;
            deltaOld[i][j] = deltaStart;
        }
    }
}

void Calculation::Delete() //TODO seg fault!!!
{
    for(int i=0; i < SIZE_Y; i++)
    {
        delete deltaNew[i];
        delete deltaOld[i];
    }
    delete deltaNew;
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
    //TODO think about what happens if an exception is thrown
}


void Calculation::SetUpModel()
{
    if(!checkInit)
    {
        throw runtime_error("Calculation was not Initialised. Run Calculation::Init() before setting up Model\n");
    }
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
    modelSetUp = true;

    //Construct model
    model.construct();
}


complex<double> Calculation::FuncDelta(Index to, Index from)
{
//    cout << "delta1" << endl;
    int from_x = from.at(0);
    int from_y = from.at(1);
    int from_s = from.at(2);

//    cout << deltaNew[from_x][from_y] << ", "<< from_s << endl;
    switch(from_s)
    {
    case 0:
        return conj(deltaOld[from_x][from_y]);
    case 1:
        return -conj(deltaOld[from_x][from_y]);
    case 2:
        return -deltaOld[from_x][from_y];
    case 3:
        return deltaOld[from_x][from_y];
    default:
        throw runtime_error("Something went wrong in Calculation::FuncDelta");
    }

    return deltaNew[from_x][from_y]*2.0*(0.5-from_s);
//    int to_x = to.at(0);
//    int to_y = to.at(1);
//    int to_s = to.at(2);
}



void Calculation::ScLoop(bool writeEachDelta)
{
//    #pragma omp parallel for
//    complex<double> pairFunction = pe.calculateExpectationValue({x,y,3},{x,y,0});

//    delta_new[x][y] = -interactionValueBCS*pairFunction;

    if(!checkInit  & !modelSetUp)
    {
        throw runtime_error("Calculation was not Initialised or model is not set up. Run Calculation::Init() and Calculation::SetUpModel() before running Calculation::ScLoop\n");
    }
    if(!cSolver)
    {
        cSolver = unique_ptr<ChebyshevSolver>( new ChebyshevSolver());
        cSolver->setModel(&model);
        cSolver->setScaleFactor(SCALE_FACTOR);
    }


    if(!pe)
    {
        pe = unique_ptr<CPropertyExtractor>(new CPropertyExtractor(cSolver.get(), //TODO check what CPropertyExtractor does with the pointer
                    NUM_COEFFICIENTS,
                    ENERGY_RESOLUTION,
                    false,
                    false,
                    true,
                    -SCALE_FACTOR,
                    SCALE_FACTOR));
    }


    int loopCounter = 0;
    complex<double>** deltaTmp;
    while(loopCounter < numberSCRuns)
    {



        deltaTmp = deltaOld;
        deltaOld = deltaNew;
        deltaNew = deltaTmp;
        double diffDelta = 0;
        double diffDeltaMax = 0;
        for(int x=0; x < SIZE_X; x++)
        {
            for(int y=0; y < SIZE_Y; y++)
            {
                deltaNew[x][y] = pe->calculateExpectationValue({x,y,3},{x,y,0})*couplingPotential;
                diffDelta = abs(deltaNew[x][y]-deltaOld[x][y]);
                if(diffDelta > diffDeltaMax)
                {
                    diffDeltaMax = diffDelta;
                }
            }
        }

        //Switching the delta in the buffer

        if(!verbose)
        {
            cout << "Loop number: " << loopCounter +1 << ", eps= " << diffDeltaMax << endl;
        }
        if(writeEachDelta)
        {
            stringstream loopFileName;
            stringstream ss;
            loopFileName << "DeltaLoop_" << loopCounter +1;
            ss << "Delta of loop number " << loopCounter +1;
            vector<complex<double>> deltaOutput = Convert2DArrayTo1DVector(deltaNew, SIZE_X, SIZE_Y);
            const int RANK = 1;
            int dims[RANK] = {SIZE_X};
            FileWriter::setFileName(loopFileName.str());
            FileWriter::clear();
            FileWriter::write(&real(deltaOutput)[0], RANK, dims, ss.str());
        }
        if(diffDeltaMax < epsDelta)
        {
            break;
            if(!verbose)
            {
                cout << "End of SC loop after " << loopCounter + 1 <<
                " iterations, delta DeltaMax: " << diffDeltaMax << endl;
            }
        }
        loopCounter++;
    }
}


vector<complex<double>> Calculation::Convert2DArrayTo1DVector(complex<double>** const input, int sizeX, int sizeY)
{
    vector<complex<double>> out;
    for(int i=0; i < sizeX; i++)
    {
        for(int j=0; j < sizeY; j++)
        {
            out.push_back(input[i][j]);
        }
    }
    return out;
}

complex<double>** Calculation::Convert1DVectorTo2DArray(vector<complex<double>> const input, int sizeX, int sizeY) //TODO needs more testing
{
    complex<double>** out;
    out = new complex<double>*[sizeX];

    for(int i=0; i < sizeX; i++)
    {
        out[i] = new complex<double>[sizeY];
        for(int j=0; j < sizeY; j++)
        {
            out[i][j] = input[j+i*sizeY];
        }
    }
    return out;
}


void Calculation::CalcLDOS()
{
    if(!checkInit  & !modelSetUp)
    {
        throw runtime_error("Calculation was not Initialised or model is not set up. Run Calculation::Init() and Calculation::SetUpModel() before running Calculation::CalcLDOS\n");
    }
    if(!cSolver)
    {
        cSolver = unique_ptr<ChebyshevSolver>( new ChebyshevSolver());
        cSolver->setModel(&model);
        cSolver->setScaleFactor(SCALE_FACTOR);
    }


    if(!pe)
    {
        //Create PropertyExtractor. The parameter are in order: The
        //ChebyshevSolver, number of expansion coefficients used in the
        //Cebyshev expansion, energy resolution with which the Green's function
        // is evaluated, whether calculate expansion functions using a GPU or
        //not, whether to evaluate the Green's function using a GPU or not,
        //whether to use a lookup table for the Green's function or not
        //(required if the Green's function is evaluated on a GPU), and the
        //lower and upper bound between which the Green's function is evaluated
        //(has to be inside the interval [-SCALE_FACTOR, SCALE_FACTOR]).
        pe = unique_ptr<CPropertyExtractor>(new CPropertyExtractor(cSolver.get(), //TODO check what CPropertyExtractor does with the pointer
                    NUM_COEFFICIENTS,
                    ENERGY_RESOLUTION,
                    false,
                    false,
                    true,
                    -SCALE_FACTOR,
                    SCALE_FACTOR));
    }

    //Extract local density of states and write to file
    Property::LDOS *ldos = pe->calculateLDOS({IDX_X, SIZE_Y/2, IDX_SUM_ALL},
                        {SIZE_X, 1, 4});



    //Set filename and remove any file already in the folder
    FileWriter::setFileName(fileName);
    FileWriter::clear();
    FileWriter::writeLDOS(ldos);
    delete ldos;
}
