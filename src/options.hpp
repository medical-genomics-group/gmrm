#pragma once
#include <string>
#include <vector>

class Options {

public:
    void read_command_line_options(int argc, char** argv);
    void print_phenotypes(void);

private:
    std::vector<std::string> phenotypes;

};
