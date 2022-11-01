#ifndef result_H
#define result_H

#include <cmath>
#include <string>
#include <vector>
using namespace std;

class Result {
    public:
        Result(int M, int k, double beta, double gamma, int no_test, double mean, double variance);
        string to_line();
        bool test_combine_res(Result r);
        Result combine_res(Result r);
    private:
        int M, k, no_test;
        double beta, gamma, mean, variance;
};

class ResultEnsemble {
    public:
        void add_result(Result r);
        void add_result(int M, int k, double beta, double gamma, int no_test, double mean, double variance);
        void to_file(const string file_name);
    private:
        vector<Result> results;
};

#endif