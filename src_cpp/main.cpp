#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <numeric>
#include <map>
#include <fstream>

#include "circle.h"
#include "grid.h"
#include "result.h"


string generate_id() {
    string chars("0123456789abcdefghijklmnopqrstuvwxyz");
    string id(12, ' ');
    for (int i = 0; i < 12; i++)
	    id[i] = chars[rand() % chars.size()];
    return id;
}

vector<vector<Point>> configuration_model(int M, int deg_v, int deg_c, double gamma) {
    vector<Point> qbts;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            qbts.push_back(Point(i, j));
        }
    }
    int N = pow(M, 2);
    int num_checks = int(N*deg_v)/deg_c;
    double r = sqrt(2)*pow(M/2, gamma);
    int c_ind, v_ind;

    vector<int> c_inds, v_inds;
    vector<int> vs(N, deg_v);
    mat<bool>* pot_qbts = new mat<bool>(num_checks, N, true);
    vector<vector<Point>> ops(num_checks, vector<Point>());


    while (accumulate(vs.begin(), vs.end(), 0)) {
        if (pot_qbts->count_nonzero()) {
            c_inds = (*pot_qbts).where();
            c_ind = c_inds[rand() % c_inds.size()];
        } else {
            cout << "Failed" << endl;
            break;
        }

        v_inds = (*pot_qbts).where(c_ind);

        v_ind = v_inds[rand() % v_inds.size()];


        ops[c_ind].push_back(qbts[v_ind]);

        // cout << "here" << endl;


        if (int(ops[c_ind].size()) == deg_c) {
            for (int i = 0; i < N; i++) {
                (*pot_qbts)(c_ind, i) = false;
            }
        } else {
            (*pot_qbts)(c_ind, v_ind) = false;
        }

        for (int i = 0; i < N; i++) {
            if ((*pot_qbts)(c_ind, i)) {
                vector<Point> tmp_op(ops[c_ind]);
                tmp_op.push_back(qbts[i]);
                Circle circ = Circle(tmp_op);
                if (circ.r > r) (*pot_qbts)(c_ind, i) = false;
            }
        }

        vs[v_ind]--;
        if (!vs[v_ind]) {
            for (int i = 0; i < num_checks; i++) {
                (*pot_qbts)(i, v_ind) = false;
            }
        }
    }

    for (size_t i = 0; i < ops.size(); i++) {
        sort(begin(ops[i]), end(ops[i]));
    }

    return ops;
}

vector<vector<Point>> randomly_draw_generators(int M, int k, double beta, double gamma) {
    vector<vector<Point>> gens;
    vector<Point> all_points;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            all_points.push_back(Point(i, j));
        }
    }

    int N = int(pow(M, 2*beta));
    double L = sqrt(2)*pow(M/2, gamma);

    while (int(gens.size()) < N) {
        Point point = all_points[rand() % all_points.size()];
        int cx = point.x, cy = point.y;
        vector<Point> in_points;
        for (size_t i = 0; i < all_points.size(); i++) {
            if (Circle::is_in_circle(Circle(cx, cy, L), all_points[i])) in_points.push_back(all_points[i]);
        }

        if (int(in_points.size()) >= k) {
            vector<Point> selected_points;
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
    if (argc != 3) return 1;
    int M = stoi(argv[1]);
    int k = stoi(argv[2]);
    int no_test = 1;
    
    vector<double> betas;
    for (int i = 10; i < 11; i++) {
        betas.push_back(i*0.1);
    }
    vector<double> gammas;
    for (int i = 10; i < 11; i++) {
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
                // vector<vector<Point>> gens = randomly_draw_generators(M, k, beta, gamma);
                vector<vector<Point>> gens = configuration_model(M, 5, 5, 1);


                // map<Point, int> counts;

                // for (size_t i = 0; i < gens.size(); i++) {
                //     for (Point qbt : gens[i]) {
                //         counts[qbt]++;
                //     }
                // }

                // for (const auto& [qbt, count] : counts) {
                //     cout << "(" << qbt.x << " " << qbt.y << "): " << count << endl;
                // }
                // for (size_t i = 0; i < gens.size(); i++) {
                //     vector<Point> gen = gens[i];
                //     for (auto it = begin(gen); it != end(gen); it++) {
                //         cout << "(" << (*it).x << " " << (*it).y << ")";
                //     }
                //     cout << endl;
                // }
                // cout << endl;

                // int rounds = grid.greedy_route_set(gens);
                int rounds = grid.route_independent_sets(gens);

                cout << rounds << endl;
                // if (rounds == -1) {
                //     cout << "0" << endl;
                //     cout << beta << " " << gamma << endl; 
                //     for (size_t i = 0; i < gens.size(); i++) {
                //         vector<Point> gen = gens[i];
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