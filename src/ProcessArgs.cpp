#include "ProcessArgs.h"

#include <iostream>
#include <exception>


using namespace std;

const int ProcessArgs::EXPECT_INPUT_FILE_PATH = -1;
const string ProcessArgs::FLAG_RESTART = "-r";
const string ProcessArgs::FLAG_INPUT_FILE_PATH = "-i";


ProcessArgs::ProcessArgs(int argc, char **argv)
{
    bool expectStringArg = false;
    int expectArgument = 0;


    for(int i=0; i < argc; i++)
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
            default:
            {
                exit(-1); //TODO error msg
                break;
            }
            }
            expectStringArg = false;
            expectArgument = 0;
        }
        else if(input.compare(FLAG_RESTART))
        {
            restart = true;
        }
        else if(input.compare(FLAG_INPUT_FILE_PATH))
        {
            expectStringArg = true;
            expectArgument = EXPECT_INPUT_FILE_PATH;
        }
        else
        {
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
