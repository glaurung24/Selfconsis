#ifndef CALCULATION_H
#define CALCULATION_H

#include <complex>
#include "Model.h"

using namespace std;
using namespace TBTK;

const complex<double> i(0, 1);


class Calculation
{
    public:
        Calculation();
        Calculation(std::string);
        virtual ~Calculation();
        void SetUpModel();
    protected:
    private:
        int N = 20;
        int SIZE_X = 4*N;
        int SIZE_Y = 2*N+1;
        int SPIN_D = 4;

        complex<double> mu = -4.0;
        complex<double> t = 1.0;
        complex<double> z = 0.5; //Zeeman coupling
        complex<double>*** delta;
        complex<double> deltaStart = 0.3;
        complex<double> alpha=0.3;

        Model model;

        int NUM_COEFFICIENTS = 1000;
        int ENERGY_RESOLUTION = 2000;
        double SCALE_FACTOR = 10.;
        complex<double>*** delta_current;
        complex<double>*** delta_old;
        complex<double>*** delta_tmp;

        complex<double> funcDelta(Index, Index);
        void InitDelta();
};

#endif // CALCULATION_H
