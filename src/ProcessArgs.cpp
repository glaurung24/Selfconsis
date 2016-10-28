#include "ProcessArgs.h"

#include <iostream>
#include <exception>


using namespace std;

const int ProcessArgs::EXPECT_INPUT_FILE_PATH = -1;
const int ProcessArgs::EXPECT_OUTPUT_FILE_PATH = -2;
const string ProcessArgs::FLAG_RESTART = "-r";
const string ProcessArgs::FLAG_INPUT_FILE_PATH = "-i";
const string ProcessArgs::FLAG_OUTPUT_FILE_PATH = "-o";


ProcessArgs::ProcessArgs(int argc, char **argv)
{
    bool expectStringArg = false;
    int expectArgument = 0;


    for(int i=1; i < argc; i++)
    {
        string input(argv[i]);
        if(expectStringArg)
        {
            switch(expectArgument)
            {
            case EXPECT_INPUT_FILE_PATH:
            {
                inputFilePath = argv[i];
                break;
            }
            case EXPECT_OUTPUT_FILE_PATH:
            {
                outputFilePath = argv[i];
                break;
            }
            default:
            {
                cout << "error in ProcessArgs" << endl;

                exit(-1); //TODO error msg
                break;
            }
            }
            expectStringArg = false;
            expectArgument = 0;
        }
        else if(input.compare(FLAG_RESTART) == 0)
        {
            restart = true;
        }
        else if(input.compare(FLAG_INPUT_FILE_PATH) == 0)
        {
            expectStringArg = true;
            expectArgument = EXPECT_INPUT_FILE_PATH;
        }
        else if(input.compare(FLAG_OUTPUT_FILE_PATH) == 0)
        {
            expectStringArg = true;
            expectArgument = EXPECT_OUTPUT_FILE_PATH;
        }
        else
        {
            cout << "error in ProcessArgs" << endl;
            exit(-1); //TODO error msg
        }
    }

}

ProcessArgs::~ProcessArgs()
{
    //dtor
}

string ProcessArgs::getInputFilePath() const
{
    return inputFilePath;
}

bool ProcessArgs::getRestart() const
{
    return restart;
}
