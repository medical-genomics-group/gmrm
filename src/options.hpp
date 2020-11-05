#pragma once
#include <string>
#include <vector>

class Options {

public:
    Options() = default;
    Options(int argc, char** argv) {
        read_command_line_options(argc, argv);
        check_options();
    }
    void read_command_line_options(int argc, char** argv);
    std::string get_bed_file() const { return bed_file; }
    std::string get_dim_file() const { return dim_file; }
    const std::vector<std::string>& get_phen_files() const { return phen_files; }
    void list_phen_files() const;
    int  count_phen_files() const { return phen_files.size(); }
    int  get_verbosity() const { return verbosity; };
    bool verbosity_level(const int level) const { return level > get_verbosity() ? false : true; }

private:
    std::string bed_file = "";
    std::string dim_file = "";
    int verbosity = 0;
    std::vector<std::string> phen_files;
    void check_options();
    void fail_if_last(char** argv, const int i);
};
