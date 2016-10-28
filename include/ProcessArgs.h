#ifndef PROCESSARGS_H
#define PROCESSARGS_H

#include <vector>
#include <tuple>
#include <string>

class ProcessArgs
{
    public:
        ProcessArgs(int argc, char **argv);
        virtual ~ProcessArgs();
        std::string getInputFilePath() const;
        std::string getOutputFilePath() const;
        bool getRestart() const;
    protected:
    private:
        static const int EXPECT_INPUT_FILE_PATH;
        static const int EXPECT_OUTPUT_FILE_PATH;
        static const std::string FLAG_RESTART;
        static const std::string FLAG_INPUT_FILE_PATH;
        static const std::string FLAG_OUTPUT_FILE_PATH;
        bool restart = false;
        std::string inputFilePath = "input";
        std::string outputFilePath = "TBTKResults.h5";
};

#endif // PROCESSARGS_H
