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
Matrix<complex<double>> Calculation::deltaNew;
Matrix<complex<double>> Calculation::deltaOld;
bool Calculation::checkInit = false;
bool Calculation::modelSetUp = false;
int Calculation::numberSCRuns;
double Calculation::epsDelta;
bool Calculation::verbose;
bool Calculation::largerBorders = false;
bool Calculation::useChebyChev = true;
bool Calculation::useGPU;
int Calculation::sCLoopCounter = 0;
bool Calculation::sCLoop = true;
unique_ptr<Util::ParameterSet> Calculation::ps;

unique_ptr<ChebyshevSolver> Calculation::cSolver;
unique_ptr<CPropertyExtractor> Calculation::pe;

string Calculation::inputFileName;
string Calculation::outputFileName;


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
const string Calculation::INIT_DELTA_REAL_ID = "DeltaReal";
const string Calculation::DELTA_LOOP_REAL_ID = "DeltaLoopReal";
const string Calculation::INIT_DELTA_IMAG_ID = "DeltaImag";
const string Calculation::DELTA_LOOP_IMAG_ID = "DeltaLoopImag";
const string Calculation::EPS_DELTA_ID = "EpsDelta";
const string Calculation::IS_MAGNETIZED_ID = "IsMagnetized";
const string Calculation::LARGER_BORDERS_ID = "DoubleBorders";
const string Calculation::USE_CHEBYCHEV_ID = "UseChebyChev";

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
    InitIsMagnetized(true);

    epsDelta = 0.05;
    numberSCRuns = 2;

    cSolver = nullptr;
    pe = nullptr;

    NUM_COEFFICIENTS = 1000;
    ENERGY_RESOLUTION = 2000;
    SCALE_FACTOR = 10.;
    sCLoopCounter = 0;

    outputFileName = "TBTKResults.h5";
    useGPU = false;


    FileWriter::setFileName(outputFileName);
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

void Calculation::InitIsMagnetized(bool magnetized)
{
    isMagnetized.reserve(SIZE_X);


    for(int i =0; i < SIZE_X; i++)
    {
        vector<bool> row(SIZE_Y, false);
        isMagnetized.push_back( row );
    }
    cout << SIZE_X << ", " << SIZE_Y << endl;

    if(magnetized)
    {
        if(largerBorders)
        {
            for(int x = SIZE_X/3; x < 2*SIZE_X/3; x++)
            {
                isMagnetized[x][SIZE_Y/2]=true;
            }
        }
        else
        {
            for(int x = SIZE_X/4; x < 3*SIZE_X/4; x++)
            {
                isMagnetized[x][SIZE_Y/2]=true;
            }
        }
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
    inputFileName = input_file;
    ps = unique_ptr<Util::ParameterSet>(FileParser::readParameterSet(inputFileName));
     //Zeeman coupling
//    counter_z = ps->getInt("counter_z");
    checkInit = true;
    if(ps->boolExists(LARGER_BORDERS_ID))
    {
        largerBorders = ps->getBool(LARGER_BORDERS_ID);
    }
    N = ps->getInt(SIZE_N_ID);
    if(largerBorders)
    {
        SIZE_X = 6*N;
        SIZE_Y = 4*N+1;
    }
    else
    {
        SIZE_X = 4*N;
        SIZE_Y = 2*N+1;
    }
    SPIN_D = 4;

    mu = ps->getComplex(CHEM_POT_ID);
    t = ps->getComplex(HOPPING_POT_ID);
    z = ps->getComplex(ZEEMAN_POT_ID); //Zeeman coupling 0.5

    deltaStart = ps->getComplex(DELTA_START_ID);
    alpha = ps->getComplex(RASHBA_COUPLING_ID);
    couplingPotential = ps->getComplex(COUPLING_POT_ID);
    periodicBoundCond = ps->getBool(PERIODIC_BOUND_ID);

    InitDelta();
    InitIsMagnetized(ps->getBool(IS_MAGNETIZED_ID));



    epsDelta = ps->getDouble(EPSILON_DELTA_ID);
    numberSCRuns = ps->getInt(MAX_NR_SCL_RUNS_ID);

    cSolver = nullptr;
    pe = nullptr;

    NUM_COEFFICIENTS = ps->getInt(NR_CHEBYCHEV_COEFF_ID);
    ENERGY_RESOLUTION = ps->getInt(ENERGY_RESOLUTION_ID);
    SCALE_FACTOR = ps->getDouble(SCALE_FACTOR_ID);
    sCLoop = ps->getBool(SC_LOOP_ID);
    sCLoopCounter = 0;
    if(ps->intExists(SC_LOOP_NR_ID))
    {
        if(ps->getInt(SC_LOOP_NR_ID))
        {
            cout << "Calculation was already started once. Use -r flag to restart." << endl;
            exit(-1);
        }
    }
    else
    {
        ps->addInt(SC_LOOP_NR_ID, sCLoopCounter);
    }


    outputFileName = ps->getString(OUTPUT_FILE_PATH_ID);
    useGPU = ps->getBool(USE_GPU_ID);

    //TODO add if scloop here...
    FileWriter::setFileName(outputFileName);
    FileReader::setFileName(outputFileName);
    FileWriter::clear();
    FileWriter::writeParameterSet(ps.get());
}

void Calculation::InitRestart(string input_file)
{
    inputFileName = input_file;
    unique_ptr<Util::ParameterSet> psInput = unique_ptr<Util::ParameterSet>(FileParser::readParameterSet(inputFileName));
    outputFileName = psInput->getString(OUTPUT_FILE_PATH_ID);
    FileWriter::setFileName(outputFileName);
    FileReader::setFileName(outputFileName);
    ps = unique_ptr<Util::ParameterSet>(FileReader::readParameterSet());
    checkInit = true;
    if(ps->boolExists(LARGER_BORDERS_ID))
    {
        largerBorders = ps->getBool(LARGER_BORDERS_ID);
    }
    N = ps->getInt(SIZE_N_ID);
    if(largerBorders)
    {
        SIZE_X = 6*N;
        SIZE_Y = 4*N+1;
    }
    else
    {
        SIZE_X = 4*N;
        SIZE_Y = 2*N+1;
    }
    SPIN_D = 4;

    mu = ps->getComplex(CHEM_POT_ID);
    t = ps->getComplex(HOPPING_POT_ID);
    z = ps->getComplex(ZEEMAN_POT_ID); //Zeeman coupling 0.5

    deltaStart = ps->getComplex(DELTA_START_ID);
    alpha = ps->getComplex(RASHBA_COUPLING_ID);
    couplingPotential = ps->getComplex(COUPLING_POT_ID);
    periodicBoundCond = ps->getBool(PERIODIC_BOUND_ID);

    InitIsMagnetized(ps->getBool(IS_MAGNETIZED_ID));

    epsDelta = ps->getDouble(EPSILON_DELTA_ID);
    numberSCRuns = ps->getInt(MAX_NR_SCL_RUNS_ID);

    cSolver = nullptr;
    pe = nullptr;

    NUM_COEFFICIENTS = ps->getInt(NR_CHEBYCHEV_COEFF_ID);
    ENERGY_RESOLUTION = ps->getInt(ENERGY_RESOLUTION_ID);
    SCALE_FACTOR = ps->getDouble(SCALE_FACTOR_ID);
    sCLoop = psInput->getBool(SC_LOOP_ID);
    sCLoopCounter = psInput->getInt(SC_LOOP_NR_ID);

    if(sCLoop)
    {
        if(!ps->intExists(SC_LOOP_NR_ID))
        {
            ps->addInt(SC_LOOP_NR_ID, sCLoopCounter);
        }
        if(!sCLoopCounter)
        {
            cout << "The SC loop count has to be at least one for restarting the calculation." << endl;
            exit(-1);
        }
        InitDelta(sCLoopCounter);
    }

    useGPU = ps->getBool(USE_GPU_ID);

    readDelta(sCLoopCounter);
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

    if(writeEachDelta && !sCLoopCounter)
    {
        WriteDelta(sCLoopCounter, -1);
        //TODO write sc loop number to file
    }

    if(!sCLoop)
    {
        if(verbose)
        {
            cout << "Omiting SC loop..." << endl;
        }
        return;
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

        double deltaRel = RelDiffDelta();
        if(verbose)
        {
            cout << "Loop number: " << sCLoopCounter +1 << ", rel eps= " << deltaRel << endl;
        }
        if(writeEachDelta)
        {
            WriteDelta(sCLoopCounter + 1, deltaRel);
            ps->setInt(SC_LOOP_NR_ID, sCLoopCounter + 1);
            FileParser::writeParameterSet(ps.get(), inputFileName);
        }

        if(deltaRel < epsDelta)
        {

            if(verbose)
            {
                cout << "End of SC loop after " << sCLoopCounter + 1 <<
                " iterations, reltive delta Delta: " << deltaRel << endl;
            }
            ps->setBool(SC_LOOP_ID, false);
            FileParser::writeParameterSet(ps.get(), inputFileName);
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

void Calculation::WriteDelta(int loopNr, double epsDelta)
{
    stringstream loopFileNameReal;
    stringstream loopFileNameImag;
//    stringstream loopFileNameArg;
    if(loopNr < 0)
    {
        loopFileNameReal << INIT_DELTA_REAL_ID;
        loopFileNameImag << INIT_DELTA_IMAG_ID;
    }
    else
    {
        if(loopNr < 10)
        {
            loopFileNameReal << DELTA_LOOP_REAL_ID << "_0" << loopNr;
            loopFileNameImag << DELTA_LOOP_IMAG_ID << "_0" << loopNr;
        }
        else
        {
            loopFileNameReal << DELTA_LOOP_REAL_ID <<  "_" << loopNr;
            loopFileNameImag << DELTA_LOOP_IMAG_ID <<  "_" << loopNr;
        }
    }

    vector<complex<double>> deltaOutput = ConvertMatrixToVector(deltaNew);


    const int RANK = 2;
    int dims[RANK] = {SIZE_X, SIZE_Y};
    FileWriter::write(&GetRealVec(deltaOutput)[0], RANK, dims, loopFileNameReal.str());
    FileWriter::write(&GetImagVec(deltaOutput)[0], RANK, dims, loopFileNameImag.str());
//    FileWriter::writeAttributes(&epsDelta, &EPS_DELTA_ID, 1, loopFileNameReal.str());
//    FileWriter::writeAttributes(&epsDelta, &EPS_DELTA_ID, 1, loopFileNameImag.str());
}

void Calculation::readDelta(int nr_sc_loop)
{
    stringstream loopFileNameReal;
    stringstream loopFileNameImag;
//    stringstream loopFileNameArg;
    if(nr_sc_loop < 10)
    {
        loopFileNameReal << DELTA_LOOP_REAL_ID << "_0" << nr_sc_loop;
        loopFileNameImag << DELTA_LOOP_IMAG_ID << "_0" << nr_sc_loop;
    }
    else
    {
        loopFileNameReal << DELTA_LOOP_REAL_ID <<  "_" << nr_sc_loop;
        loopFileNameImag << DELTA_LOOP_IMAG_ID <<  "_" << nr_sc_loop;
    }


    double* delta_real_from_file = nullptr;
    double* delta_imag_from_file = nullptr;
    int rank;
    int* dims;
    FileReader::read(&delta_real_from_file, &rank, &dims, loopFileNameReal.str());
    delete [] dims;
    FileReader::read(&delta_imag_from_file, &rank, &dims, loopFileNameImag.str());

    Matrix<double> realPart = ConvertVectorToMatrix(delta_real_from_file, SIZE_X, SIZE_Y);
    Matrix<double> imagPart = ConvertVectorToMatrix(delta_imag_from_file, SIZE_X, SIZE_Y);

    for(unsigned int i=0; i < deltaNew.size(); i++)
    {
        for(unsigned int j=0; j < deltaNew[i].size(); j++)
        {
            deltaNew[i][j] = realPart[i][j] + i*imagPart[i][j];//TODO complex constructor calling with argument and absolute value
        }
    }
    delete [] dims;
    delete [] delta_real_from_file;
    delete [] delta_imag_from_file;
}

vector<double> Calculation::GetRealVec(vector<complex<double>> input)
{
    vector<double> output;
    output.reserve(input.size());
    for(unsigned int i=0; i < input.size(); i++)
    {
        output.push_back(real(input[i]));
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
            diffDelta = abs(real(deltaNew[x][y]-deltaOld[x][y]));
            if(diffDelta > diffDeltaMax)
            {
                diffDeltaMax = diffDelta;
                deltaRel = diffDeltaMax/abs(real(deltaOld[x][y]));
            }
        }
    }
    return deltaRel;
}

vector<double> Calculation::GetImagVec(vector<complex<double>> input)
{
    vector<double> output;
    for(unsigned int i=0; i < input.size(); i++)
    {
        output.push_back(imag(input[i]));
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

Matrix<double> Calculation::ConvertVectorToMatrix(const double *input, int sizeX, int sizeY)
{
    Matrix<double> out;
    out.reserve(sizeX);
    for(int i=0; i < sizeX; i++)
    {
        vector<double> row;
        row.reserve(sizeY);
        for(int j=0; j < sizeY; j++)
        {
            row.push_back(input[j+i*sizeY]);
        }
        out.push_back(row);
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
                        {SIZE_X, 1, 4});

    //Set filename and remove any file already in the folder

    FileWriter::writeLDOS(ldos, "LDOS");
    delete ldos;
}

void Calculation::setVerbose(bool input)
{
    verbose = input;
}

void Calculation::setLargerBorders(bool input)
{
    largerBorders = input;
}

