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
    protected:
    private:
        static const int EXPECT_INPUT_FILE_PATH;
        static const std::string FLAG_RESTART;
        static const std::string FLAG_INPUT_FILE_PATH;
        bool restart = false;
        std::string inputFilePath = "input";


        std::vector<std::tuple<std::string, bool>> flags;

};

#endif // PROCESSARGS_H
