#ifndef CALCULATION_H
#define CALCULATION_H

#include <complex>
#include <memory>
#include "Model.h"
#include "ChebyshevSolver.h"
#include "CPropertyExtractor.h"

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


    protected:
    private:
        Calculation();
        virtual ~Calculation();


        static void InitDelta();
        static complex<double> FuncDelta(Index, Index);
        static vector<complex<double>> Convert2DArrayTo1DVector(complex<double>** const, int, int);
        static complex<double>** Convert1DVectorTo2DArray(vector<complex<double>> const, int, int);
        static vector<double> getAbsVec(vector<complex<double>>);
        static vector<double> getPhaseVec(vector<complex<double>>);
        static void setBoundary(complex<double>**);



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

        static Model model;

        static int NUM_COEFFICIENTS;
        static int ENERGY_RESOLUTION;
        static double SCALE_FACTOR;
        static complex<double>** deltaNew;
        static complex<double>** deltaOld;
        static bool checkInit;
        static bool modelSetUp;

        static int numberSCRuns;
        static int epsDelta;
        static bool verbose;
        static unique_ptr<ChebyshevSolver> cSolver;
        static unique_ptr<CPropertyExtractor> pe;
        static string fileName;
};

#endif // CALCULATION_H
