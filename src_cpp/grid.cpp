#include "grid.h"

Grid::Grid(int inp_M, int inp_k) {
    M = inp_M;
    k = inp_k;

    ancillas = new mat<int>(M, M, k);
    dests = new mat<bool>(M, M, false);
    in_progress = new mat<int>(M, M, 0);

    vector<Generator> generators;
    vector<Generator> finished_gens;
    vector<vector<vector<Point>>> full_chains;
    vector<vector<tuple<Point, Point>>> bell_pairs;
}

Grid::~Grid() {
    delete ancillas;
    delete dests;
    delete in_progress;
}

const mat<int>& Grid::get_ancillas() {
    return *ancillas;
}

vector<Point> Grid::find_chain(const mat<int>& grid, Point site1, Point site2) {
    mat<bool> visited(M, M, false);
    mat<Point> parent(M, M, Point());
    deque<Point> queue = {site1};
    visited(site1.x, site1.y) = true;

    if ((grid(site1.x, site1.y) <= 0) || 
        (grid(site2.x, site2.y) <= 0)) {
        return {};
    }

    while (queue.size()) {
        Point curr_site = queue.front();
        int x = curr_site.x, y = curr_site.y;
        queue.pop_front();

        if (curr_site == site2) {
            vector<Point> chain = {curr_site};
            while (curr_site != site1) {
                curr_site = parent(curr_site.x, curr_site.y);
                chain.insert(begin(chain), curr_site);
            }
            return chain;
        }

        if ((curr_site != site1) && (grid(x, y) < 2)) { 
            continue; 
        }
        vector<Point> pot_nbrs = { Point(x, y+1), Point(x, y-1), Point(x+1, y), Point(x-1, y) };
        for (auto it = begin(pot_nbrs); it != end(pot_nbrs); it++) {
            int new_x = (*it).x, new_y = (*it).y;

            if ((0 <= new_x && new_x < M) && (0 <= new_y && new_y < M)) {
                if (!visited(new_x, new_y)) {
                    parent(new_x, new_y) = curr_site;
                    queue.push_back(*it);
                    visited(new_x, new_y) = true;
                }
            }
        }

    }

    return {};
}

void Grid::add_chain(vector<Point> chain, int gen) {
    for (size_t i = 0; i < chain.size(); i++) {
        int x = chain[i].x, y = chain[i].y;
        if ((i == 0) || (i == (chain.size() - 1))) {
            (*ancillas)(x, y) = (*ancillas)(x, y) - 1;
        } else {
            (*ancillas)(x, y) = (*ancillas)(x, y) - 2;
        }
    }
    full_chains[gen].push_back(chain);
}

void Grid::tmp_add_chain(mat<int>& grid, vector<Point> chain) {
    for (size_t i = 0; i < chain.size(); i++) {
        int x = chain[i].x, y = chain[i].y;
        if ((i == 0) || (i == (chain.size() - 1))) {
            grid(x, y) = grid(x, y) - 1;
        } else {
            grid(x, y) = grid(x, y) - 2;
        }
    }
}

void Grid::perform_bell_measurement() {
    for (size_t i = 0; i < generators.size(); i++) {
        vector<vector<Point>> gen_chains = full_chains[generators[i].key];
        for (size_t j = 0; j < gen_chains.size(); j++) {
            vector<Point> chain = gen_chains[j];
            for (size_t k = 1; k < (chain.size() - 1); k++) {
                int x = chain[k].x, y = chain[k].y;
                (*ancillas)(x, y) = (*ancillas)(x, y) + 2;
            }
            bell_pairs[generators[i].key].push_back(make_tuple(chain[0], chain[chain.size()-1]));
        } 
        full_chains[generators[i].key] = vector<vector<Point>>();
    }
}

void Grid::perform_syndrome_measurements() {
    for (size_t i = 0; i < generators.size(); i++) {
        Generator gen = generators[i];
        if (gen.is_done()) {
            int x = gen.dest.x, y = gen.dest.y;
            (*ancillas)(x, y) = (*ancillas)(x, y) + 1;
            (*dests)(x, y) = false;

            vector<tuple<Point, Point>> pairs = bell_pairs[gen.key];
            for (size_t j = 0; j < pairs.size(); j++) {
                const auto[qbt1, qbt2] = pairs[j];
                int x1 = qbt1.x, y1 = qbt1.y;
                (*ancillas)(x1, y1) = (*ancillas)(x1, y1) + 1;
                int x2 = qbt2.x, y2 = qbt2.y;
                (*ancillas)(x2, y2) = (*ancillas)(x2, y2) + 1;
            }

            bell_pairs[gen.key] = vector<tuple<Point, Point>>();

            vector<Point> qbts = gen.all_qbts;
            for (size_t j = 0; j < qbts.size(); j++) {
                Point qbt = qbts[j];
                (*in_progress)(qbt.x, qbt.y) = (*in_progress)(qbt.x, qbt.y) - 1;
            }
        }
    }
}

void Grid::reset_generator(Generator& gen) {
    if ((!gen.is_done()) && (gen.dest != Point())) { // not done but started routing
        int x = gen.dest.x, y = gen.dest.y;
        (*ancillas)(x, y) = (*ancillas)(x, y) + 1;
        (*dests)(x, y) = false;

        vector<tuple<Point, Point>> pairs = bell_pairs[gen.key];
        for (size_t j = 0; j < pairs.size(); j++) {
            const auto[qbt1, qbt2] = pairs[j];
            int x1 = qbt1.x, y1 = qbt1.y;
            (*ancillas)(x1, y1) = (*ancillas)(x1, y1) + 1;
            int x2 = qbt2.x, y2 = qbt2.y;
            (*ancillas)(x2, y2) = (*ancillas)(x2, y2) + 1;
        }

        bell_pairs[gen.key] = vector<tuple<Point, Point>>();

        vector<Point> qbts = gen.all_qbts;
        for (size_t j = 0; j < qbts.size(); j++) {
            Point qbt = qbts[j];
            (*in_progress)(qbt.x, qbt.y) = (*in_progress)(qbt.x, qbt.y) - 1;
        }
    }

    gen.reset();
}

vector<vector<Point>> Grid::route_generator(vector<Point> gen, Point prior_dest) {
    bool is_prior_dest = (prior_dest != Point());
    vector<vector<vector<Point>>> out = {};
    vector<Point> possible_dests = gen;
    if (is_prior_dest) {
        possible_dests = {prior_dest};
    }

    for (size_t i = 0; i < possible_dests.size(); i++) {
        Point dest = possible_dests[i];
        int x = dest.x, y = dest.y;

        // prevent gridlock, could have parameter turning this off
        if (!is_prior_dest && ((*in_progress)(x, y))) continue;

        mat<int> tmp_grid = get_ancillas();
        vector<vector<Point>> chains = {};
        vector<Point> routed_qbts = {};
        size_t tot_len = 0;

        // dest is the meeting site
        if (!is_prior_dest) {
            tmp_grid(x, y) = tmp_grid(x, y) - 1;
        }

        for (size_t j = 0; j < gen.size(); j++) {
            Point site = gen[j];
            
            if ((dest != site) && (!(*dests)(site.x, site.y))) {
                vector<Point> chain = find_chain(tmp_grid, dest, site);
                if (chain.size()) {
                    tmp_add_chain(tmp_grid, chain);
                    chains.push_back(chain);
                    routed_qbts.push_back(site);
                    tot_len = tot_len + chain.size();
                } else continue;
            }
        }
        out.push_back(chains);
    }

     if (!is_prior_dest) {
        stable_sort(begin(out), end(out), []
            (const vector<vector<Point>>& lhs, const vector<vector<Point>>& rhs) {
                if (lhs.size() < rhs.size()) return false;
                if (rhs.size() < lhs.size()) return true;

                size_t tot_lhs_size = 0;
                size_t tot_rhs_size = 0;
                for (size_t i = 0; i < lhs.size(); i++) {
                    tot_lhs_size += lhs[i].size();
                }
                for (size_t i = 0; i < rhs.size(); i++) {
                    tot_rhs_size += rhs[i].size();
                }
                if (tot_lhs_size < tot_rhs_size) return true;
                if (tot_rhs_size < tot_lhs_size) return false;
                return false;
            });
    }

    // for (size_t i = 0; i < out.size(); i++) {
    //     vector<vector<Point>> chains = out[i];
    //     for (size_t j = 0; j < chains.size(); j++) {
    //         vector<Point> chain = chains[j];
    //         for (auto it = begin(chain); it != end(chain); it++) {
    //             cout << "(" << (*it).x << " " << (*it).y << ")";
    //         }
    //         cout << endl;
    //     }
    //     cout << "........." << endl;
    // }

    if (out.size()) {
        return out[0];
    } else {
        return {};
    }
}

int Grid::greedy_route_set(vector<vector<Point>> gens) {
    for (size_t i = 0; i < gens.size(); i++) {
        if (int(gens[i].size()) > k) return 0;

        full_chains.push_back(vector<vector<Point>>());
        bell_pairs.push_back(vector<tuple<Point, Point>>());

        Generator gen(gens[i], i);
        generators.push_back(gen);
    }

    // stable_sort(begin(generators), end(generators), []
    //     (const Generator& lhs, const Generator& rhs) {
    //         return lhs.c.r < rhs.c.r;
    //     });

    for (size_t rounds = 0; rounds < (gens.size() + 1); rounds++) {
        stable_sort(begin(generators), end(generators), []
            (const Generator& lhs, const Generator& rhs) {
                if (lhs.num_routed() == rhs.num_routed()) {
                    // return lhs.c.r < rhs.c.r;
                    return false;
                } else {
                    return lhs.num_routed() > rhs.num_routed();
                }
            });

        for (size_t i = 0; i < generators.size(); i++) {
            Generator &gen = generators[i];
            vector<vector<Point>> chains = route_generator(gen.qbts_to_route, gen.dest);

            if (chains.size()) {
                if (gen.dest == Point()) {
                    Point dest = chains[0][0];
                    int x = dest.x, y = dest.y;
                    gen.dest = dest;
                    (*ancillas)(x, y) = (*ancillas)(x, y) - 1;
                    (*dests)(x, y) = true;
                    gen.route_qbt(dest);
                    gen.start = rounds;

                    vector<Point> qbts = gen.all_qbts;
                    for (size_t j = 0; j < qbts.size(); j++) {
                        Point qbt = qbts[j];
                        (*in_progress)(qbt.x, qbt.y) = (*in_progress)(qbt.x, qbt.y) + 1;
                    }
                }

                for (size_t j = 0; j < chains.size(); j++) {
                    vector<Point> chain = chains[j];
                    gen.route_qbt(chain[chain.size() - 1]);
                    add_chain(chain, gen.key);
                }
            }      
            // reset_generator(gen);

        }

        // double sum = 0;
        // for (int i = 0; i < M; i++) {
        //     for (int j = 0; j < M; j++) {
        //         sum += (*ancillas)(i, j);
        //     }
        // }
        // cout << rounds << "," << sum << "," << M*M*k << "," << sum/double(M*M*k) << endl;
          
        perform_bell_measurement();
        perform_syndrome_measurements();

        for (size_t i = 0; i < generators.size(); i++) {
            Generator gen = generators[i];
            if (gen.is_done()) {
                cout << gen.c.r << "," << rounds+1 << "," << rounds+1-gen.start << endl;   
            }
        }

        generators.erase(remove_if(begin(generators), end(generators), 
            [](const Generator& gen) {
                return gen.is_done();
            }), end(generators));
        if (!generators.size()) {
            return rounds+1;
        }        
    }

    return -1;
}


void Grid::greedy_route_reccuring(vector<vector<Point>> gens, int max_rounds, int restart) {
    for (size_t i = 0; i < gens.size(); i++) {
        if (int(gens[i].size()) > k) return;

        full_chains.push_back(vector<vector<Point>>());
        bell_pairs.push_back(vector<tuple<Point, Point>>());

        Generator gen(gens[i], i);
        generators.push_back(gen);
    }

    // what info do we want out of this function? 
    // every time a generator is restarted, print round, cycle, and if it finished...
    // maybe when adding back in you should shuffle the generators
    for (int rounds = 1; rounds < max_rounds; rounds++) {
        stable_sort(begin(generators), end(generators), []
            (const Generator& lhs, const Generator& rhs) {
                if (lhs.num_routed() == rhs.num_routed()) {
                    // return lhs.c.r < rhs.c.r;
                    return false;
                } else {
                    return lhs.num_routed() > rhs.num_routed();
                }
            });

        for (size_t i = 0; i < generators.size(); i++) {
            Generator &gen = generators[i];
            vector<vector<Point>> chains = route_generator(gen.qbts_to_route, gen.dest);

            if (chains.size()) {
                if (gen.dest == Point()) {
                    Point dest = chains[0][0];
                    int x = dest.x, y = dest.y;
                    gen.dest = dest;
                    (*ancillas)(x, y) = (*ancillas)(x, y) - 1;
                    (*dests)(x, y) = true;
                    gen.route_qbt(dest);

                    vector<Point> qbts = gen.all_qbts;
                    for (size_t j = 0; j < qbts.size(); j++) {
                        Point qbt = qbts[j];
                        (*in_progress)(qbt.x, qbt.y) = (*in_progress)(qbt.x, qbt.y) + 1;
                    }
                }

                for (size_t j = 0; j < chains.size(); j++) {
                    vector<Point> chain = chains[j];
                    gen.route_qbt(chain[chain.size() - 1]);
                    add_chain(chain, gen.key);
                }
            }      
        }

        

        perform_bell_measurement();
        perform_syndrome_measurements();

        double sum = 0;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < M; j++) {
                sum += (*ancillas)(i, j);
            }
        }


        for (auto it = begin(generators); it != end(generators);) {
            Generator &gen = *it;

            if (gen.is_done()) {
                finished_gens.push_back(gen);
                it = generators.erase(it);
            } else {
                it++;
            }
        }

        if ((rounds) && (rounds % restart == 0)) {
            cout << "ROUNDS:" << rounds << "," << sum << "," << M*M*k << "," << sum/double(M*M*k) << endl;

            for (auto it = begin(generators); it != end(generators); it++) {
                Generator &gen = *it;
                cout << gen.c.r << "," << gen.key << "," << rounds << "," << gen.cycle << "," << gen.is_done() << endl;
            }
            for (auto it = begin(finished_gens); it != end(finished_gens);) {
                Generator &gen = *it;

                cout << gen.c.r << "," << gen.key << "," << rounds << "," << gen.cycle << "," << gen.is_done() << endl;
                reset_generator(gen);

                generators.push_back(gen);
                it = finished_gens.erase(it);
            }
        }
    }
}


int Grid::route_independent_sets(vector<vector<Point>> gens) {
    // split the generators up into indepenent sets to be routed separately.
    // each independent set can be routed in one round. Total time is thus number of independent sets

    size_t orig_size = gens.size();
    for (size_t rounds = 0; rounds < (orig_size + 1); rounds++) {
        vector<vector<Point>> ind_set = {gens[0]};
        gens.erase(begin(gens));

        for (auto it = begin(gens); it != end(gens);) {
            bool add = true;

            for (auto it2 = begin(ind_set); it2 != end(ind_set); it2++) {
                vector<Point> intersection;
                set_intersection(begin(*it), end(*it),
                               begin(*it2), end(*it2),
                               back_inserter(intersection));

                if (intersection.size()) {
                    add = false;
                }
            }
            
            if (add) {
                ind_set.push_back(*it);
                it = gens.erase(it);
            } else {
                it++;
            }
        }

        if (!gens.size()) {
            return rounds+1;
        }
    }

    return -1;
}