#include <iostream>
#include <mpi.h>
#include <limits.h>
#include "utilities.hpp"


double round_dp(const double in) {
    return in;

    printf("in = %20.15f\n", in);
    double out = round(in * 1.0E12) / 1.0E12;
    printf("ou = %20.15f\n", out);
    return out;
}

void check_malloc(const void* ptr, const int linenumber, const char* filename) {
    if (ptr == NULL) {
        fprintf(stderr, "#FATAL#: malloc failed on line %d of %s\n", linenumber, filename);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

void check_mpi(const int error, const int linenumber, const char* filename) {
    if (error != MPI_SUCCESS) {
        fprintf(stderr, "*FATAL*: MPI error %d at line %d of file %s\n", error, linenumber, filename);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

// Check whether a size_t can be casted to int or would overflow
int check_int_overflow(const size_t n, const int linenumber, const char* filename) {

    if (n > INT_MAX) {
        fprintf(stderr, "FATAL  : integer overflow detected on line %d of %s. %lu does not fit in type int.\n", linenumber, filename, n);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    return static_cast<int>(n);
}

/**
* @autor Diego Garcia
* Split a string in n parts given a splitting character.
* @param line: the string to split.
* @param charSplit: the splitting character (e.g. tab '\t', space ' ').
* @return a vector of strings with the parts of the splitted string.
*/
std::vector<std::string> split_string(std::string line, char charSplit) {
    std::vector<std::string> tokens;
    std::istringstream iss(line);
    std::string token;
    while (std::getline(iss, token, charSplit))   // but we can specify a different one
        tokens.push_back(token);

    return tokens;
}

/**
* @autor Diego Garcia
* Counts the number of lines of a file.
* @param file_path: a string with the full path of a file.
* @return an integer with the number of lines of the file.
*/
int get_file_line_count(std::string file_path) {
    int count = 0;
    std::string s;
    std::ifstream in;
    in.open(file_path);

    while (!in.eof()) {
        getline(in, s);
        if (s.length() == 0)  continue;
        count++;
    }

    in.close();
    return count;
}

/**
* @autor Diego Garcia
* Given a effects matrix whose values are -1, 0, or 1, this function
* converts those values to 0, 1, or 2, respectively
* @param matrix_effects: a vector of vectors of strings with the effect matrix.
* @return a vector of vectors of integers with the parsed effect matrix.
*/
std::vector <std::vector<int>> parse_effect_values(std::vector <std::vector<std::string>> matrix_effects) {
    std::vector <std::vector<int>> matrix_effects_ardyh(matrix_effects.size());
    std::vector<std::string> tmp_line;
    std::vector<int> tmp_line_ardyh;

    for (int i = 0; i < matrix_effects.size(); i++) {
        tmp_line = matrix_effects[i];
        tmp_line_ardyh = std::vector<int>(tmp_line.size());
        for (int j = 0; j < tmp_line.size(); j++) {
            if (tmp_line[j] == "-1") { tmp_line_ardyh[j] = 0; } else
            if (tmp_line[j] == "0") { tmp_line_ardyh[j] = 1; }  else
            if (tmp_line[j] == "1") { tmp_line_ardyh[j] = 2; }
        }

        matrix_effects_ardyh[i] = tmp_line_ardyh;
    }

    return matrix_effects_ardyh;
}
