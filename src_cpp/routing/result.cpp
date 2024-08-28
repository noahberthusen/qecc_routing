#include "result.h"
#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>


Result::Result(int M, int k, double beta, double gamma, int no_test, double mean, double variance) :
    M(M), k(k), no_test(no_test), beta(beta), gamma(gamma), mean(mean), variance(variance) {}

bool Result::test_combine_res(Result r) {
    return (M == r.M && k == r.k && 
        (abs(beta - r.beta) < 0.00001) && (abs(gamma - r.gamma) < 0.00001));
}

Result Result::combine_res(Result r) {
    // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm
    // https://online.stat.psu.edu/stat414/lesson/24/24.4
    int new_no_test = no_test + r.no_test;
    double new_mean = mean + ((r.mean - mean) / new_no_test);
    double new_variance = variance + ((r.mean - mean) * (r.mean - new_mean));
    return Result(M, k, beta, gamma, new_no_test, new_mean, new_variance);
}

string Result::to_line() {
    stringstream oss;
    oss << M << "," << k << "," << beta << "," << gamma << "," << no_test << "," << mean << "," << variance/(no_test-1) << "\n";
    return oss.str();
}

void ResultEnsemble::add_result(Result r) {
    for (vector<Result>::iterator it=begin(results); it != end(results); it++) {
        if (it->test_combine_res(r)) {
            *it = it->combine_res(r);
            return;
        }
    }
    results.push_back(r);
}

void ResultEnsemble::add_result(int M, int k, double beta, double gamma, int no_test, double mean, double variance) {
    Result r(M, k, beta, gamma, no_test, mean, variance);
    add_result(r);
}

void ResultEnsemble::to_file(const string file_name) {
    vector<string> lines(results.size());
    for (size_t i = 0; i < results.size(); i++) {
        lines[i] = results[i].to_line();
    }
    sort(begin(lines), end(lines));

    const string tmp_file_name = file_name + ".tmp";
    ofstream file(tmp_file_name, ios::out | ios::trunc);

    file << "m,k,beta,gamma,no_test,mean,variance\n";
    for (vector<string>::iterator it=begin(lines); it != end(lines); it++) {
        file << *it;
    }
    file.close();
    rename(tmp_file_name.c_str(), file_name.c_str());
}