#ifndef CALCULATION_H
#define CALCULATION_H

#include <complex>
#include <memory>
#include "Model.h"
#include "ChebyshevSolver.h"
#include "CPropertyExtractor.h"
#include "Matrix.h"
#include "ParameterSet.h"

using namespace std;
using namespace TBTK;

const complex<double> i(0, 1);


class Calculation
{
    public:
        static void Init();
        static void Init(string);
        static void InitRestart(string);
        static void SetUpModel();
        static void ScLoop(bool writeEachDelta = false);
        static void Delete();
        static void CalcLDOS();
        static void setVerbose(bool);
        static void setLargerBorders(bool);


    protected:
    private:
        Calculation();
        virtual ~Calculation();


        static void InitDelta();
        static void InitDelta(int);
        static complex<double> FuncDelta(Index, Index);
        static vector<complex<double>> ConvertMatrixToVector(const Matrix<complex<double>>&);
        static Matrix<double> ConvertVectorToMatrix(const double *, int , int );
//        static complex<double>** Convert1DVectorTo2DArray(vector<complex<double>> const, int, int);
        static vector<double> GetRealVec(vector<complex<double>>);
        static vector<double> GetImagVec(vector<complex<double>>);
        static void SetBoundary(complex<double>**);
        static double RelDiffDelta();
        static void WriteDelta(int, double);
        static void SwapDeltas();
        static void InitIsMagnetized(bool);
        static void readDelta(int);



        static int N;
        static int SIZE_X;
        static int SIZE_Y;
        static int SPIN_D;

        static complex<double> mu;
        static complex<double> t;
        static complex<double> z; //Zeeman coupling
        static complex<double> deltaStart;
        static complex<double> alpha;
        static complex<double> couplingPotential;
        static Matrix<bool> isMagnetized;
        static bool periodicBoundCond;

        static Model model;

        static int NUM_COEFFICIENTS;
        static int ENERGY_RESOLUTION;
        static double SCALE_FACTOR;
        static Matrix<complex<double>> deltaNew;
        static Matrix<complex<double>>  deltaOld;
        static bool checkInit;
        static bool modelSetUp;
        static bool sCLoop;
        static int sCLoopCounter;
        static unique_ptr<ParameterSet> ps;

        static int numberSCRuns;
        static double epsDelta;
        static bool verbose;
        static bool largerBorders;
        static bool useChebyChev;
        static unique_ptr<ChebyshevSolver> cSolver;
        static unique_ptr<CPropertyExtractor> pe;
        static string inputFileName;
        static string outputFileName;
        static bool useGPU;

        static const string SIZE_N_ID;
        static const string CHEM_POT_ID;
        static const string HOPPING_POT_ID;
        static const string ZEEMAN_POT_ID;
        static const string DELTA_START_ID;
        static const string RASHBA_COUPLING_ID;
        static const string COUPLING_POT_ID;
        static const string PERIODIC_BOUND_ID;
        static const string EPSILON_DELTA_ID;
        static const string MAX_NR_SCL_RUNS_ID;
        static const string NR_CHEBYCHEV_COEFF_ID;
        static const string ENERGY_RESOLUTION_ID;
        static const string SCALE_FACTOR_ID;
        static const string USE_GPU_ID;
        static const string OUTPUT_FILE_PATH_ID;
        static const string SC_LOOP_ID;
        static const string SC_LOOP_NR_ID;
        static const string INIT_DELTA_REAL_ID;
        static const string DELTA_LOOP_REAL_ID;
        static const string INIT_DELTA_IMAG_ID;
        static const string DELTA_LOOP_IMAG_ID;
        static const string EPS_DELTA_ID;
        static const string IS_MAGNETIZED_ID;
        static const string LARGER_BORDERS_ID;
        static const string USE_CHEBYCHEV_ID;
};

#endif // CALCULATION_H
