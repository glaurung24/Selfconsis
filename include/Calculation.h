#ifndef CALCULATION_H
#define CALCULATION_H

#include <complex>
#include <memory>
#include "Model.h"
#include "ChebyshevSolver.h"
#include "CPropertyExtractor.h"
#include "Matrix.h"

using namespace std;
using namespace TBTK;

const complex<double> i(0, 1);


class Calculation
{
    public:
        static void Init();
        static void Init(string);
        static void SetUpModel();
        static void ScLoop(bool writeEachDelta = false);
        static void Delete();
        static void CalcLDOS();
        static void SetZeemanPot(complex<double>);


    protected:
    private:
        Calculation();
        virtual ~Calculation();


        static void InitDelta();
        static complex<double> FuncDelta(Index, Index);
        static vector<complex<double>> ConvertMatrixToVector(const Matrix<complex<double>>&);
//        static complex<double>** Convert1DVectorTo2DArray(vector<complex<double>> const, int, int);
        static vector<double> GetAbsVec(vector<complex<double>>);
        static vector<double> GetPhaseVec(vector<complex<double>>);
        static void SetBoundary(complex<double>**);
        static double RelDiffDelta();
        static void WriteDelta(int);
        static void SwapDeltas();
        static void InitIsMagnetized();



        static int N;
        static int SIZE_X;
        static int SIZE_Y;
        static int SPIN_D;

        static complex<double> mu;
        static complex<double> t;
        static complex<double> z; //Zeeman coupling
        static complex<double>*** delta;
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

        static int numberSCRuns;
        static double epsDelta;
        static bool verbose;
        static unique_ptr<ChebyshevSolver> cSolver;
        static unique_ptr<CPropertyExtractor> pe;
        static string fileName;
        static bool useGPU;
};

#endif // CALCULATION_H
