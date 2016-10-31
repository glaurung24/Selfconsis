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
#include "FileParser.h"
#include "FileReader.h"


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
Matrix<bool> Calculation::isMagnetized;
bool Calculation::periodicBoundCond;

Model Calculation::model;

int Calculation::NUM_COEFFICIENTS;
int Calculation::ENERGY_RESOLUTION;
double Calculation::SCALE_FACTOR;
Matrix<complex<double>>  Calculation::deltaNew;
Matrix<complex<double>> Calculation::deltaOld;
bool Calculation::checkInit = false;
bool Calculation::modelSetUp = false;
int Calculation::numberSCRuns;
double Calculation::epsDelta;
bool Calculation::verbose;
bool Calculation::useGPU;
int Calculation::sCLoopCounter = 0;
bool Calculation::sCLoop = true;
unique_ptr<Util::ParameterSet> Calculation::ps;

unique_ptr<ChebyshevSolver> Calculation::cSolver;
unique_ptr<CPropertyExtractor> Calculation::pe;

string Calculation::fileName;


const string Calculation::SIZE_N_ID = "SizeN";
const string Calculation::CHEM_POT_ID = "ChemPot";
const string Calculation::HOPPING_POT_ID = "HoppingPot";
const string Calculation::ZEEMAN_POT_ID = "ZeemanPot";
const string Calculation::DELTA_START_ID = "DeltaStart";
const string Calculation::RASHBA_COUPLING_ID = "RashbaCoupling";
const string Calculation::COUPLING_POT_ID = "CouplingPot";
const string Calculation::PERIODIC_BOUND_ID = "PeriodicBound";
const string Calculation::EPSILON_DELTA_ID = "EpsDelta";
const string Calculation::MAX_NR_SCL_RUNS_ID = "MaxNrSCRuns";
const string Calculation::NR_CHEBYCHEV_COEFF_ID = "NumChebCoef";
const string Calculation::ENERGY_RESOLUTION_ID = "EnergyRes";
const string Calculation::SCALE_FACTOR_ID = "ScaleFactor";
const string Calculation::USE_GPU_ID = "UseGPU";
const string Calculation::SC_LOOP_ID = "SCLoop";
const string Calculation::SC_LOOP_NR_ID = "SCLoopNr";
const string Calculation::OUTPUT_FILE_PATH_ID = "OutputFilePath";
const string Calculation::INIT_DELTA_ABS_ID = "DeltaAbs";
const string Calculation::DELTA_LOOP_ABS_ID = "DeltaLoopAbs";
const string Calculation::INIT_DELTA_ARG_ID = "DeltaArg";
const string Calculation::DELTA_LOOP_ARG_ID = "DeltaLoopArg";




void Calculation::Init()
{
    checkInit = true;
    N = 2;
    SIZE_X = 4*N;
    SIZE_Y = 2*N+1;
    SPIN_D = 4;

    mu = -4.0;
    t = 1.0;
    z = 0.5; //Zeeman coupling 0.5

    deltaStart = 0.3;
    alpha = 0.3;
    couplingPotential = 5;
    periodicBoundCond = false;

    InitDelta();
    InitIsMagnetized();

    epsDelta = 0.05;
    numberSCRuns = 2;

    cSolver = nullptr;
    pe = nullptr;

    NUM_COEFFICIENTS = 1000;
    ENERGY_RESOLUTION = 2000;
    SCALE_FACTOR = 10.;
    sCLoopCounter = 0;

    fileName = "TBTKResults.h5";
    useGPU = false;


    FileWriter::setFileName(fileName);
    FileWriter::clear();

}

void Calculation::InitDelta()
{
    deltaNew.reserve(SIZE_X);
    deltaOld.reserve(SIZE_X);


    for(int i =0; i < SIZE_X; i++)
    {
        vector<complex<double>> row(SIZE_Y, deltaStart);
        deltaNew.push_back( row );
        deltaOld.push_back( row );
//        deltaNew[i].reserve();
//        deltaOld[i].reserve(SIZE_Y);
//        for(unsigned int j=0; j < deltaNew[i].size(); j++)
//        {
//            deltaNew[i].push_back(deltaStart); //TODO mabye faster way?
//            deltaOld[i].push_back(deltaStart);
//        }
    }
}

void Calculation::InitDelta(int nr_sc_loop)
{
    deltaNew.reserve(SIZE_X);
    deltaOld.reserve(SIZE_X);


    for(int i =0; i < SIZE_X; i++)
    {
        vector<complex<double>> row(SIZE_Y, deltaStart);
        deltaNew.push_back( row );
        deltaOld.push_back( row );
    }
    readDelta(nr_sc_loop);
}

void Calculation::InitIsMagnetized()
{
    isMagnetized.reserve(SIZE_X);


    for(int i =0; i < SIZE_X; i++)
    {
        vector<bool> row(SIZE_Y, false);
        isMagnetized.push_back( row );
    }

    for(int x = SIZE_X/4; x < 3*SIZE_X/4; x++){
        isMagnetized[x][SIZE_Y/2]=true;
    }
}


void Calculation::Delete()
{
//    for(int i=0; i < SIZE_Y; i++)
//    {
//        delete deltaNew[i];
//        delete deltaOld[i];
//    }
//    delete deltaNew;
//    delete deltaOld;
}



void Calculation::Init(std::string input_file) //TODO
{

    ps = unique_ptr<Util::ParameterSet>(FileParser::readParameterSet(input_file));
     //Zeeman coupling
//    counter_z = ps->getInt("counter_z");
    checkInit = true;
    N = ps->getInt(SIZE_N_ID);
    SIZE_X = 4*N;
    SIZE_Y = 2*N+1;
    SPIN_D = 4;

    mu = ps->getComplex(CHEM_POT_ID);
    t = ps->getComplex(CHEM_POT_ID);
    z = ps->getComplex(ZEEMAN_POT_ID); //Zeeman coupling 0.5

    deltaStart = ps->getComplex(DELTA_START_ID);
    alpha = ps->getComplex(RASHBA_COUPLING_ID);
    couplingPotential = ps->getComplex(COUPLING_POT_ID);
    periodicBoundCond = ps->getBool(PERIODIC_BOUND_ID);

    InitDelta();
    InitIsMagnetized();

    epsDelta = ps->getDouble(EPSILON_DELTA_ID);
    numberSCRuns = ps->getInt(MAX_NR_SCL_RUNS_ID);

    cSolver = nullptr;
    pe = nullptr;

    NUM_COEFFICIENTS = ps->getInt(NR_CHEBYCHEV_COEFF_ID);
    ENERGY_RESOLUTION = ps->getInt(ENERGY_RESOLUTION_ID);
    SCALE_FACTOR = ps->getDouble(SCALE_FACTOR_ID);
    sCLoop = ps->getBool(SC_LOOP_ID);
    sCLoopCounter = 0;

    if(sCLoop)
    {
        ps->addInt(SC_LOOP_NR_ID, sCLoopCounter);
    }

    fileName = ps->getString(OUTPUT_FILE_PATH_ID);
    useGPU = ps->getBool(USE_GPU_ID);

    //TODO add if scloop here...


    FileWriter::setFileName(fileName);
    FileReader::setFileName(fileName);
    FileWriter::clear();
    FileWriter::writeParameterSet(ps.get());
}



void Calculation::InitRestart(string output_file)
{

    ps = unique_ptr<Util::ParameterSet>(FileReader::readParameterSet(output_file));
    checkInit = true;
    N = ps->getInt(SIZE_N_ID);
    SIZE_X = 4*N;
    SIZE_Y = 2*N+1;
    SPIN_D = 4;

    mu = ps->getComplex(CHEM_POT_ID);
    t = ps->getComplex(CHEM_POT_ID);
    z = ps->getComplex(ZEEMAN_POT_ID); //Zeeman coupling 0.5

    deltaStart = ps->getComplex(DELTA_START_ID);
    alpha = ps->getComplex(RASHBA_COUPLING_ID);
    couplingPotential = ps->getComplex(COUPLING_POT_ID);
    periodicBoundCond = ps->getBool(PERIODIC_BOUND_ID);

    InitIsMagnetized();

    epsDelta = ps->getDouble(EPSILON_DELTA_ID);
    numberSCRuns = ps->getInt(MAX_NR_SCL_RUNS_ID);

    cSolver = nullptr;
    pe = nullptr;

    NUM_COEFFICIENTS = ps->getInt(NR_CHEBYCHEV_COEFF_ID);
    ENERGY_RESOLUTION = ps->getInt(ENERGY_RESOLUTION_ID);
    SCALE_FACTOR = ps->getDouble(SCALE_FACTOR_ID);
    sCLoop = ps->getBool(SC_LOOP_ID);

    if(sCLoop)
    {
        sCLoopCounter = ps->getInt(SC_LOOP_NR_ID);
        if(!sCLoopCounter)
        {
            cout << "The SC loop count has to be at least one for restarting the calculation." << endl;
            exit(-1);
        }
        InitDelta(sCLoopCounter);
    }

    fileName = ps->getString(OUTPUT_FILE_PATH_ID);
    useGPU = ps->getBool(USE_GPU_ID);


    //TODO set delta


    //TODO do some more stuff here

    FileWriter::setFileName(output_file);
    FileReader::setFileName(fileName);


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
                if(periodicBoundCond || x+1 < SIZE_X){
                    model.addHAAndHC(HoppingAmplitude(-t,	{(x+1)%SIZE_X, y, s},	{x, y, s}));
                    model.addHAAndHC(HoppingAmplitude(t,	{x, y, s+2},{(x+1)%SIZE_X, y, s+2}));
//					model.addHAAndHC(HoppingAmplitude(t, {(x+1)%SIZE_X, y, s+2}, {x, y, s+2})); //same results with this line as with the line above
                }
                if(periodicBoundCond || y+1 < SIZE_Y){
                    model.addHAAndHC(HoppingAmplitude(-t,	{x, (y+1)%SIZE_Y, s},	{x, y, s}));
                    model.addHAAndHC(HoppingAmplitude(t,  {x, y, s+2}, {x, (y+1)%SIZE_Y, s+2}));
//					model.addHAAndHC(HoppingAmplitude(t,   {x, (y+1)%SIZE_Y, s+2},{x, y, s+2}));
                }

//------------------------Rashba hopping term--------------------------------------

                if(periodicBoundCond || x+1 < SIZE_X){
    //                    model.addHAAndHC(HoppingAmplitude(alpha*2.0*(0.5-s), {(x+1)%SIZE_X, y, s*2},	{x, y, s*2+1}));
    //                    model.addHAAndHC(HoppingAmplitude(-alpha*2.0*(0.5-s), {x, y, s*2},	{(x+1)%SIZE_X, y, s*2+1}));
                    model.addHAAndHC(HoppingAmplitude(alpha *2.0*(0.5-s), {(x+1)%SIZE_X,y,(s+1)%2}, {x,y,s}));
                    model.addHAAndHC(HoppingAmplitude(-alpha *2.0*(0.5-s), {x,y,s+2}, {(x+1)%SIZE_X,y,(s+1)%2+2}));
                }

                if(periodicBoundCond || y+1 < SIZE_Y){
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
    if(useGPU)
    {
        model.constructCOO();
    }
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
                    useGPU,
                    false,
                    true,
                    -SCALE_FACTOR,
                    SCALE_FACTOR));
    }

    //TODO add if scloop here...




    if(writeEachDelta)
    {
        WriteDelta(sCLoopCounter);
    }

    while(sCLoopCounter < numberSCRuns)
    {



        SwapDeltas();
        if(useGPU)
        {
            model.reconstructCOO();
        }


        #pragma omp parallel for
        for(int x=0; x < SIZE_X; x++)
        {
            for(int y=0; y < SIZE_Y; y++)
            {
                deltaNew[x][y] = -pe->calculateExpectationValue({x,y,3},{x,y,0})*couplingPotential;

            }
        }


//        setBoundary(deltaNew);

        double deltaRel = RelDiffDelta();
        if(!verbose)
        {
            cout << "Loop number: " << sCLoopCounter +1 << ", rel eps= " << deltaRel << endl;
        }
        if(writeEachDelta)
        {
            WriteDelta(sCLoopCounter + 1);
        }

        if(deltaRel < epsDelta)
        {

            if(!verbose)
            {
                cout << "End of SC loop after " << sCLoopCounter + 1 <<
                " iterations, reltive delta Delta: " << deltaRel << endl;
            }
            //TODO write to paramset SC=false
            break;
        }
        sCLoopCounter++;
    }
}

void Calculation::SwapDeltas()
{
        Matrix<complex<double>> deltaTmp;
        deltaTmp = deltaOld;
        deltaOld = deltaNew;
        deltaNew = deltaTmp;
}

void Calculation::WriteDelta(int loopNr)
{
    stringstream loopFileNameAbs;
    stringstream loopFileNameArg;
//    stringstream loopFileNameArg;
    if(loopNr < 0)
    {
        loopFileNameAbs << INIT_DELTA_ABS_ID;
        loopFileNameArg << INIT_DELTA_ARG_ID;
    }
    else
    {
        if(loopNr < 10)
        {
            loopFileNameAbs << DELTA_LOOP_ARG_ID << "_0" << loopNr;
            loopFileNameArg << DELTA_LOOP_ARG_ID << "_0" << loopNr;
        }
        else
        {
            loopFileNameAbs << DELTA_LOOP_ABS_ID <<  "_" << loopNr;
            loopFileNameArg << DELTA_LOOP_ARG_ID <<  "_" << loopNr;
        }
    }

    vector<complex<double>> deltaOutput = ConvertMatrixToVector(deltaNew);


    const int RANK = 2;
    int dims[RANK] = {SIZE_X, SIZE_Y};
    FileWriter::write(&GetAbsVec(deltaOutput)[0], RANK, dims, loopFileNameAbs.str());
    FileWriter::write(&GetPhaseVec(deltaOutput)[0], RANK, dims, loopFileNameArg.str());
}

void Calculation::readDelta(int nr_sc_loop)
{
    stringstream loopFileNameAbs;
    stringstream loopFileNameArg;
//    stringstream loopFileNameArg;
    if(loopNr < 10)
    {
        loopFileNameAbs << DELTA_LOOP_ARG_ID << "_0" << loopNr;
        loopFileNameArg << DELTA_LOOP_ARG_ID << "_0" << loopNr;
    }
    else
    {
        loopFileNameAbs << DELTA_LOOP_ABS_ID <<  "_" << loopNr;
        loopFileNameArg << DELTA_LOOP_ARG_ID <<  "_" << loopNr;
    }


    double** delta_abs_from_file;
    double** delta_arg_from_file;
    int rank;
    int* dims;
    FileReader::read(delta_abs_from_file, &rank, &dims, loopFileNameAbs.str());
    delete [] dims;
    FileReader::read(delta_arg_from_file, &rank, &dims, loopFileNameAbs.str());

    for(int i=0; i < deltaNew.size(); i++)
    {
        for(int j=0; j < deltaNew[i].size(); j++)
        {
            deltaNew[i][j] = //TODO complex constructor calling with argument and absolute value
        }
    }
    delete [] dims;
}


vector<double> Calculation::GetAbsVec(vector<complex<double>> input)
{
    vector<double> output;
    output.reserve(input.size());
    for(unsigned int i=0; i < input.size(); i++)
    {
        output.push_back(abs(input[i]));
    }
    return output;
}


double Calculation::RelDiffDelta()
{
    double diffDelta = 0;
    double diffDeltaMax = 0;
    double deltaRel = 1;
    for(int x=0; x < SIZE_X; x++)
    {
        for(int y=0; y < SIZE_Y; y++)
        {
            diffDelta = abs(deltaNew[x][y]-deltaOld[x][y]);
            if(diffDelta > diffDeltaMax)
            {
                diffDeltaMax = diffDelta;
                deltaRel = diffDeltaMax/abs(deltaOld[x][y]);
            }
        }
    }
    return deltaRel;
}
vector<double> Calculation::GetPhaseVec(vector<complex<double>> input)
{
    vector<double> output;
    for(unsigned int i=0; i < input.size(); i++)
    {
        output.push_back(arg(input[i]));
    }
    return output;
}

vector<complex<double>> Calculation::ConvertMatrixToVector(const Matrix<complex<double>>& input)
{
    vector<complex<double>> out;
    for(unsigned int i=0; i < input.size(); i++)
    {
        for(unsigned int j=0; j < input[i].size(); j++)
        {
            out.push_back(input[i][j]);
        }
    }
    return out;
}

void Calculation::SetBoundary(complex<double>** input)
{
    //Upper y boundary
    for(int i=0; i < SIZE_X; i++)
    {
        input[i][SIZE_Y-1] = deltaStart;
    }

    //Lower y boundary
    for(int i=0; i < SIZE_X; i++)
    {
        input[i][0]= deltaStart;
    }

    //left x boundary
    for(int i=0; i < SIZE_Y; i++)
    {
        input[0][i] = deltaStart;
    }

    //Right x boundary
    for(int i=0; i < SIZE_Y; i++)
    {
        input[SIZE_X-1][i] = deltaStart;
    }
}

//complex<double>** Calculation::Convert1DVectorTo2DArray(vector<complex<double>> const input, int sizeX, int sizeY) //TODO outdated
//    complex<double>** out;
//    out = new complex<double>*[sizeX];
//
//    for(int i=0; i < sizeX; i++)
//    {
//        out[i] = new complex<double>[sizeY];
//        for(int j=0; j < sizeY; j++)
//        {
//            out[i][j] = input[j+i*sizeY];
//        }
//    }
//    return out;
//}


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


//    if(!pe) //TODO
//    {
//        //Create PropertyExtractor. The parameter are in order: The
//        //ChebyshevSolver, number of expansion coefficients used in the
//        //Cebyshev expansion, energy resolution with which the Green's function
//        // is evaluated, whether calculate expansion functions using a GPU or
//        //not, whether to evaluate the Green's function using a GPU or not,
//        //whether to use a lookup table for the Green's function or not
//        //(required if the Green's function is evaluated on a GPU), and the
//        //lower and upper bound between which the Green's function is evaluated
//        //(has to be inside the interval [-SCALE_FACTOR, SCALE_FACTOR]).
//        pe = unique_ptr<CPropertyExtractor>(new CPropertyExtractor(cSolver.get(), //TODO check what CPropertyExtractor does with the pointer
//                    NUM_COEFFICIENTS,
//                    ENERGY_RESOLUTION,
//                    useGPU,
//                    false,
//                    true,
//                    -1,
//                    1));
//    }

    pe = unique_ptr<CPropertyExtractor>(new CPropertyExtractor(cSolver.get(), //TODO check what CPropertyExtractor does with the pointer
                    NUM_COEFFICIENTS,
                    ENERGY_RESOLUTION,
                    useGPU,
                    false,
                    true,
                    -1,
                    1));

    //Extract local density of states and write to file
    Property::LDOS *ldos = pe->calculateLDOS({IDX_X, SIZE_Y/2, IDX_SUM_ALL},
                        {SIZE_X, 1, 2});

    //Set filename and remove any file already in the folder

    FileWriter::writeLDOS(ldos);
    delete ldos;
}

void Calculation::SetZeemanPot(complex<double> newZeemanPot)
{
    z = newZeemanPot;
}
