#include <iostream>
#include <cmath>
#include <random>
#include <string>
#include <numeric>
#include <map>
#include <algorithm>
#include <chrono>
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

vector<vector<Point>> configuration_model(int M, int deg_v, int deg_c, double gamma) {
    // can't have operators with only one qubit in them. breaks the algo, and is not realistic for the model
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
    vector<vector<Point>> ops(num_checks, vector<Point>());

    vector<vector<int>> pot_qbts(num_checks, vector<int>());
    for (int i = 0; i < num_checks; i++) {
        for (int j = 0; j < N; j++) {
            pot_qbts[i].push_back(j);
        }
        c_inds.push_back(i);
    }

    while (accumulate(vs.begin(), vs.end(), 0)) {
        if (c_inds.size()) {
            c_ind = c_inds[rand() % c_inds.size()];
        } else {
            break;
        }

        v_inds = pot_qbts[c_ind];
        v_ind = v_inds[rand() % v_inds.size()];
        ops[c_ind].push_back(qbts[v_ind]);

        if (int(ops[c_ind].size()) == deg_c) {
            pot_qbts[c_ind].clear();
            c_inds.erase(remove(begin(c_inds), end(c_inds), c_ind), end(c_inds));
        } else {
            pot_qbts[c_ind].erase(remove(begin(pot_qbts[c_ind]),
                                         end(pot_qbts[c_ind]), v_ind),
                                         end(pot_qbts[c_ind]));
            if (!pot_qbts[c_ind].size())
                c_inds.erase(remove(begin(c_inds), end(c_inds), c_ind), end(c_inds));
        }

        // for (int i = 0; i < N; i++) {
        //     if ((*pot_qbts)(c_ind, i)) {
        //         vector<Point> tmp_op(ops[c_ind]);
        //         tmp_op.push_back(qbts[i]);
        //         Circle circ = Circle(tmp_op);
        //         if (circ.r > r) (*pot_qbts)(c_ind, i) = false;
        //     }
        // }

        if (ops[c_ind].size() == 1) {
            for (int i = 0; i < N; i++) {
                if (ops[c_ind][0].distance(qbts[i]) > (1+gamma)*r) {
                    pot_qbts[c_ind].erase(remove(begin(pot_qbts[c_ind]),
                                                 end(pot_qbts[c_ind]), i),
                                                 end(pot_qbts[c_ind]));
                    if (!pot_qbts[c_ind].size())
                        c_inds.erase(remove(begin(c_inds), end(c_inds), c_ind), end(c_inds));
                }
            }
        }
        vs[v_ind]--;
        if (!vs[v_ind]) {
            for (int i = 0; i < num_checks; i++) {
                pot_qbts[i].erase(remove(pot_qbts[i].begin(),
                    pot_qbts[i].end(), v_ind), pot_qbts[i].end());
                if (!pot_qbts[i].size())
                    c_inds.erase(remove(begin(c_inds), end(c_inds), i), end(c_inds));
            }
        }

    }

    for (auto it = begin(ops); it != end(ops);) {
        if ((*it).size() == 1) {
            it = ops.erase(it);
        } else {
            it++;
        }
    }

    for (auto it = begin(ops); it != end(ops); it++) {
        sort(begin(*it), end(*it));
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

vector<vector<Point>> draw_from_distribution(int M, int k, double beta) {
    // draws from an exponential distribution
    vector<vector<Point>> gens;
    vector<Point> all_points;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < M; j++) {
            all_points.push_back(Point(i, j));
        }
    }

    int a = 2;
    auto f = [a](double x) { return a*exp(-a*x); };
    auto f_inv = [a](double x) { return -log(x/a)/a; };
    random_device rd;
    mt19937 gen(rd());
    std::uniform_real_distribution<> dist(f(0), f(1));
    // std::uniform_real_distribution<> dist(0, 1);

    int N = int(pow(M, 2*beta));

    while (int(gens.size()) < N) {
        double L = sqrt(2)*pow(M/2, f_inv(dist(gen)));
        // double L = sqrt(2)*pow(M/2, dist(gen));

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

double fn(int M) {
    // https://helloacm.com/c-function-to-compute-numerical-integral-using-function-pointers/
    int a = 2;
    auto dist = [a, M](double x) {
        return (a/(1-exp(-a)))*exp(-a*x)*pow(M,x);
    };

    double c = 1.85;
    double d = 6.35;
    double n = 1000; // increase this for more accurate integral
    // double gamma = log2(r)/log2(M/2);
    double gamma = 0.5883925047211545;

    double step = gamma / n;
    double area = 0.0;
    for (int i = 0; i < n; i++) {
        area += dist((i + 0.5) * step) * step;
    }
    return c*area + d;
}


int main(int argc, char* argv[]) {
    if (argc != 3) return 1;
    int M = stoi(argv[1]);
    int k = stoi(argv[2]);
    int no_test = 1;

    // don't do beta = 1 with configuration model
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

    vector<vector<vector<Point>>> all_gens(gammas.size(), vector<vector<Point>>());
    vector<vector<Point>> gens;

    auto start = chrono::high_resolution_clock::now();

    /*
    for (int r = 0; r < no_test; r++) {
        // if (r % int(pow(M,2)/2) == 0) {
        //     for (size_t i = 0; i < gammas.size(); i++) {
        //         all_gens[i] = configuration_model(M, k, k, gammas[i]);
        //     }
        // }

        for (size_t beta_ind = 0; beta_ind < betas.size(); beta_ind++) {
            double beta = betas[beta_ind];
            for (size_t gamma_ind = 0; gamma_ind < gammas.size(); gamma_ind++) {
                double gamma = gammas[gamma_ind];
                Grid grid(M, k+1);

                if (gens.size()) gens.clear();

                // sample(begin(all_gens[gamma_ind]),
                //     end(all_gens[gamma_ind]),
                //     back_inserter(gens),
                //     int(pow(M, 2*beta)),
                //     mt19937{random_device{}()});
                gens = draw_from_distribution(M, k, beta);

                for (auto it = begin(gens); it != end(gens);) {
                    vector<Point> gen = *it;
                    Circle c = Circle(gen);

                    if (log2(c.r)/log2(M/2) > 0.5883925047211545) {
                        it = gens.erase(it);
                    } else {
                        it++;
                    }

                }
                // gens = randomly_draw_generators(M, k, beta, gamma);

                // map<Point, int> counts;

                // for (size_t i = 0; i < gens.size(); i++) {
                //     for (Point qbt : gens[i]) {
                //         counts[qbt]++;
                //     }
                // }

                // for (const auto& [qbt, count] : counts) {
                //     cout << "(" << qbt.x << " " << qbt.y << "): " << count << endl;
                // }
                // double sum = 0;
                // for (size_t i = 0; i < gens.size(); i++) {
                //     vector<Point> gen = gens[i];
                //     // Circle c = Circle(gen);
                //     // cout << c.r << endl;
                //     for (auto it = begin(gen); it != end(gen); it++) {
                //         cout << "(" << (*it).x << " " << (*it).y << ")";
                //     }
                //     cout << endl;
                // }
                // cout << "............." << endl;
                // cout << sum/gens.size() << " " << sqrt(2)*pow(M/2, gamma);

                // cout << endl;

                grid.greedy_route_reccuring(gens, 50, fn(M));
                // int rounds = grid.route_independent_sets(gens);

                // cout << rounds << endl;
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
    */

    auto m_dist = [](Point p1, Point p2) { return abs(p1.x - p2.x) + abs(p1.y - p2.y); };

    for (M = 10; M < 101; M += 10) {
        // gens = draw_from_distribution(M, k, 1);
        gens = randomly_draw_generators(M, k, 1, 1);

        int tot_edges = 0;
        for (int _ = 0; _ < no_test; _++) {
            for (size_t i = 0; i < gens.size(); i++) {
                for (int j = 1; j < k; j++) {
                    tot_edges += m_dist(gens[i][0], gens[i][j]);
                }
            }
        }
        cout << M << "," << tot_edges/no_test << endl;
    }

    // gens = draw_from_distribution(M, k, 1);

    // int num_edges;
    // for (size_t i = 0; i < gens.size(); i++) {
    //     vector<Point> gen = gens[i];
    //     Circle c = Circle(gen);

    //     num_edges = 0;
    //     for (int j = 1; j < k; j++) {
    //         num_edges += m_dist(gen[0], gen[j]);
    //     }

    //     cout << c.r << "," << num_edges << endl;

    // }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);

    // cout << duration.count() << endl;
    return 0;
}