#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <map>
#include <fstream>
#include "grid.h"
#include "result.h"


string generate_id() {
    string chars("0123456789abcdefghijklmnopqrstuvwxyz");
    string id(12, ' ');
    for (int i = 0; i < 12; i++)
	    id[i] = chars[rand() % chars.size()];
    return id;
}

bool in_circle(tuple<int, int> point, tuple<int, int, double> circle) {
    const auto[x, y] = point;
    const auto[cx, cy, r] = circle;
    return sqrt(pow(x - cx, 2) + pow(y - cy, 2)) <= r;
}

vector<vector<tuple<int, int>>> randomly_draw_generators(int M, int k, double beta, double gamma) {
    vector<vector<tuple<int, int>>> gens;
    vector<tuple<int, int>> all_points;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            all_points.push_back(make_tuple(i, j));
        }
    }

    int N = int(pow(M, 2*beta));
    double L = sqrt(2)*pow(M, gamma);

    while (int(gens.size()) < N) {
        const auto[cx, cy] = all_points[rand() % all_points.size()];
        vector<tuple<int, int>> in_points;
        for (size_t i = 0; i < all_points.size(); i++) {
            if (in_circle(all_points[i], make_tuple(cx, cy, L))) in_points.push_back(all_points[i]);
        }

        if (int(in_points.size()) >= k) {
            vector<tuple<int, int>> selected_points;
            sample(begin(in_points), 
                end(in_points), 
                back_inserter(selected_points),
                k,
                mt19937{random_device{}()});

            sort(begin(selected_points), end(selected_points));
            gens.push_back(selected_points);
        }
    }

    return gens;
}


int main(int argc, char* argv[]) {
    int M = stoi(argv[1]);
    int k = stoi(argv[2]);
    int no_test = 1;
    
    vector<double> betas;
    for (int i = 10; i < 11; i++) {
        betas.push_back(i*0.1);
    }
    vector<double> gammas;
    for (int i = 0; i < 1; i++) {
        gammas.push_back(i*0.1);
    }

    int seed = time(NULL) + 123 * 0;
    srand (seed);
    string res_file_name = string("../results/cpp_") + generate_id() + string(".res");
    ResultEnsemble res_ens;

    for (int r = 0; r < no_test; r++) {
        for (size_t beta_ind = 0; beta_ind < betas.size(); beta_ind++) {
            double beta = betas[beta_ind];
            for (size_t gamma_ind = 0; gamma_ind < gammas.size(); gamma_ind++) {
                double gamma = gammas[gamma_ind];

                Grid grid(M, k+1);
                vector<vector<tuple<int, int>>> gens = randomly_draw_generators(M, k, beta, gamma);


                map<tuple<int, int>, int> counts;

                for (size_t i = 0; i < gens.size(); i++) {
                    for (tuple<int, int> qbt : gens[i]) {
                        counts[qbt]++;
                    }
                }

                for (const auto& [qbt, count] : counts) {
                    cout << "(" << get<0>(qbt) << " " << get<1>(qbt) << "): " << count << endl;
                }
                // for (size_t i = 0; i < gens.size(); i++) {
                //     vector<tuple<int, int>> gen = gens[i];
                //     for (auto it = begin(gen); it != end(gen); it++) {
                //         cout << "(" << get<0>(*it) << " " << get<1>(*it) << ")";
                //     }
                //     cout << endl;
                // }

                // int rounds = grid.greedy_route_set(gens);
                // int rounds = grid.route_independent_sets(gens);
                // if (rounds == -1) {
                //     cout << "0" << endl;
                //     cout << beta << " " << gamma << endl; 
                //     for (size_t i = 0; i < gens.size(); i++) {
                //         vector<tuple<int, int>> gen = gens[i];
                //         for (auto it = begin(gen); it != end(gen); it++) {
                //             cout << "(" << get<0>(*it) << " " << get<1>(*it) << ")";
                //         }
                //         cout << endl;
                //     }
                // }


                // if (rounds != -1) {
                //     res_ens.add_result(M, k, beta, gamma, 1, rounds, 0);
                // }
            }
        }
        // res_ens.to_file(res_file_name);
    }

    return 0;
}