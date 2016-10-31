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
        bool getVerbose() const;
        bool getWriteDelta() const;
    protected:
    private:
        static const int EXPECT_INPUT_FILE_PATH;
        static const int EXPECT_OUTPUT_FILE_PATH;
        static const std::string FLAG_RESTART;
        static const std::string FLAG_INPUT_FILE_PATH;
        static const std::string FLAG_OUTPUT_FILE_PATH;
        static const std::string FLAG_VERBOSE;
        static const std::string FLAG_WRITE_DELTA;
        static const std::string FLAG_NOT_WRITE_DELTA;
        bool restart = false;
        bool verbose = false;
        bool writeDelta = true;
        std::string inputFilePath = "input";
        std::string outputFilePath = "TBTKResults.h5";
};

#endif // PROCESSARGS_H
