/**
 *  @file main.cpp
 *  @brief Main function for calculating ...
 *
 *  TODO
 *
 *
 *  @author Andreas Theiler
 */

#include <complex>
#include <iostream>
#include <memory>
#include "Model.h"
#include "Util.h"
#include <sstream>
#include "H5Cpp.h" //TODO
#include <iomanip>
#include <ctime>

#include "Calculation.h"
#include "ProcessArgs.h"


using namespace std;
using namespace TBTK;


int main(int argc, char **argv){
//    time_t t = time(NULL);
//    stringstream ss;
//    ss << "TBTKLog" << put_time(localtime(&t), "%d-%m-%Y_%H-%M-%S");
//    Util::Streams::openLog(ss.str());
    Util::Streams::openLog();


    ofstream logFile("TBTKLog", fstream::out | fstream::app);
    time_t t = time(NULL);
    logFile << "\n\n###" << put_time(localtime(&t), "%d-%m-%Y_%H-%M-%S") << "\n\n";
    Util::Streams::log.rdbuf(logFile.rdbuf());

    ProcessArgs args(argc, argv);

    Calculation::Init("input");
//    Calculation::SetUpModel();
//    Calculation::ScLoop(true);
//    Calculation::CalcLDOS();

    logFile.close(); //TODO
	return 0;
}


