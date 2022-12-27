#include "grid.h"

Grid::Grid(int inp_N, int inp_k) {
    N = inp_N;
    k = inp_k;

    ancillas = new mat<int>(N, N, k);
    dests = new mat<bool>(N, N, false);
    in_progress = new mat<int>(N, N, 0);

    vector<Generator> generators;
    vector<vector<vector<tuple<int, int>>>> full_chains;
    vector<vector<tuple<tuple<int, int>, tuple<int, int>>>> bell_pairs;
}

Grid::~Grid() {
    delete ancillas;
    delete dests;
    delete in_progress;
}

const mat<int>& Grid::get_ancillas() {
    return *ancillas;
}

vector<tuple<int, int>> Grid::find_chain(const mat<int>& grid, tuple<int, int> site1, tuple<int, int> site2) {
    mat<bool> visited(N, N, false);
    mat<tuple<int,int>> parent(N, N, make_tuple(-1, -1));
    deque<tuple<int, int>> queue = {site1};
    visited(get<0>(site1), get<1>(site1)) = true;

    if ((grid(get<0>(site1), get<1>(site1)) <= 0) || 
        (grid(get<0>(site2), get<1>(site2)) <= 0)) {
        return {};
    }

    while (queue.size()) {
        tuple<int, int> curr_site = queue.front();
        const auto[y, x] = curr_site;
        queue.pop_front();

        if (curr_site == site2) {
            // for (int i = 0; i < N; i++) {
            //     for (int j = 0; j < N; j++) {
            //         tuple<int, int> site = parent(i, j);
            //         cout << "(" << get<0>(site) << "," << get<1>(site) << ")";
            //     }
            //     cout << endl;
            // }
            vector<tuple<int, int>> chain = {curr_site};
            while (curr_site != site1) {
                curr_site = parent(get<0>(curr_site), get<1>(curr_site));
                chain.insert(begin(chain), curr_site);
            }
            return chain;
        }

        if ((curr_site != site1) && (grid(y, x) < 2)) { 
            continue; 
        }
        vector<tuple<int, int>> pot_nbrs = { make_tuple(y+1, x), make_tuple(y-1, x), make_tuple(y, x+1), make_tuple(y, x-1) };
        for (auto it = begin(pot_nbrs); it != end(pot_nbrs); it++) {
            const auto[new_y, new_x] = *it;

            if ((0 <= new_x && new_x < N) && (0 <= new_y && new_y < N)) {
                if (!visited(new_y, new_x)) {
                    parent(new_y, new_x) = curr_site;
                    queue.push_back(*it);
                    visited(new_y, new_x) = true;
                }
            }
        }

    }

    return {};
}

void Grid::add_chain(vector<tuple<int, int>> chain, int gen) {
    for (size_t i = 0; i < chain.size(); i++) {
        const auto[y, x] = chain[i];
        if ((i == 0) || (i == (chain.size() - 1))) {
            (*ancillas)(y, x) = (*ancillas)(y, x) - 1;
        } else {
            (*ancillas)(y, x) = (*ancillas)(y, x) - 2;
        }
    }
    full_chains[gen].push_back(chain);
}

void Grid::tmp_add_chain(mat<int>& grid, vector<tuple<int, int>> chain) {
    for (size_t i = 0; i < chain.size(); i++) {
        const auto[y, x] = chain[i];
        if ((i == 0) || (i == (chain.size() - 1))) {
            grid(y, x) = grid(y, x) - 1;
        } else {
            grid(y, x) = grid(y, x) - 2;
        }
    }
}

void Grid::perform_bell_measurement() {
    for (size_t i = 0; i < generators.size(); i++) {
        vector<vector<tuple<int, int>>> gen_chains = full_chains[generators[i].get_key()];
        for (size_t j = 0; j < gen_chains.size(); j++) {
            vector<tuple<int, int>> chain = gen_chains[j];
            for (size_t k = 1; k < (chain.size() - 1); k++) {
                const auto[y, x] = chain[k];
                (*ancillas)(y, x) = (*ancillas)(y, x) + 2;
            }
            bell_pairs[generators[i].get_key()].push_back(make_tuple(chain[0], chain[chain.size()-1]));
        } 
        full_chains[generators[i].get_key()] = vector<vector<tuple<int, int>>>();
    }
}

void Grid::perform_syndrome_measurements() {
    for (size_t i = 0; i < generators.size(); i++) {
        Generator gen = generators[i];
        if (gen.is_done()) {
            const auto[y, x] = gen.get_dest();
            (*ancillas)(y, x) = (*ancillas)(y, x) + 1;

            vector<tuple<tuple<int, int>, tuple<int, int>>> pairs = bell_pairs[gen.get_key()];
            for (size_t j = 0; j < pairs.size(); j++) {
                const auto[qbt1, qbt2] = pairs[j];
                const auto[y1, x1] = qbt1;
                (*ancillas)(y1, x1) = (*ancillas)(y1, x1) + 1;
                const auto[y2, x2] = qbt2;
                (*ancillas)(y2, x2) = (*ancillas)(y2, x2) + 1;
            }

            bell_pairs[gen.get_key()] = vector<tuple<tuple<int, int>, tuple<int, int>>>();
        }
    }
}

vector<vector<tuple<int, int>>> Grid::route_generator(vector<tuple<int, int>> gen, tuple<int, int> prior_dest) {
    bool is_prior_dest = (prior_dest != make_tuple(-1, -1));
    vector<vector<vector<tuple<int, int>>>> out = {};
    vector<tuple<int, int>> possible_dests = gen;
    if (is_prior_dest) {
        possible_dests = {prior_dest};
    }

    for (size_t i = 0; i < possible_dests.size(); i++) {
        tuple<int, int> dest = possible_dests[i];
        const auto[y, x] = dest;

        // prevent gridlock, could have parameter turning this off
        if (!is_prior_dest && ((*in_progress)(y, x))) continue;

        mat<int> tmp_grid = get_ancillas();
        vector<vector<tuple<int, int>>> chains = {};
        vector<tuple<int, int>> routed_qbts = {};
        size_t tot_len = 0;

        // dest is the meeting site
        if (!is_prior_dest) {
            tmp_grid(y, x) = tmp_grid(y, x) - 1;
        }

        for (size_t j = 0; j < gen.size(); j++) {
            tuple<int, int> site = gen[j];
            

            if ((dest != site) && (!(*dests)(get<0>(site), get<1>(site)))) {
                vector<tuple<int, int>> chain = find_chain(tmp_grid, dest, site);
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
            (const vector<vector<tuple<int, int>>>& lhs, const vector<vector<tuple<int, int>>>& rhs) {
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
    //     vector<vector<tuple<int, int>>> chains = out[i];
    //     for (size_t j = 0; j < chains.size(); j++) {
    //         vector<tuple<int, int>> chain = chains[j];
    //         for (auto it = begin(chain); it != end(chain); it++) {
    //             cout << "(" << get<0>(*it) << " " << get<1>(*it) << ")";
    //         }
    //         cout << endl;
    //     }
    //     cout << "........." << endl;
    // }
    // cout << "++++++++" << endl;


    if (out.size()) {
        return out[0];
    } else {
        return {};
    }
}

int Grid::greedy_route_set(vector<vector<tuple<int, int>>> gens) {
    for (size_t i = 0; i < gens.size(); i++) {
        if (int(gens[i].size()) > k) return 0;

        full_chains.push_back(vector<vector<tuple<int, int>>>());
        bell_pairs.push_back(vector<tuple<tuple<int, int>, tuple<int, int>>>());

        Generator gen(gens[i], i);
        generators.push_back(gen);
    }

    for (size_t rounds = 0; rounds < (gens.size() + 1); rounds++) {
        stable_sort(begin(generators), end(generators), []
            (const Generator& lhs, const Generator& rhs) {
                if (lhs.num_routed() < rhs.num_routed()) return false;
                if (rhs.num_routed() < lhs.num_routed()) return true;
                return false;
            });

        for (size_t i = 0; i < generators.size(); i++) {
            Generator &gen = generators[i];
            vector<vector<tuple<int, int>>> chains = route_generator(gen.get_qbts_to_route(), gen.get_dest());

            if (chains.size()) {
                if (gen.get_dest() == make_tuple(-1, -1)) {
                    tuple<int, int> dest = chains[0][0];
                    const auto[y, x] = dest;
                    gen.set_dest(dest);
                    (*ancillas)(y, x) = (*ancillas)(y, x) - 1;
                    (*dests)(y, x) = true;
                    gen.route_qbt(dest);

                    vector<tuple<int, int>> qbts = gen.get_qbts();
                    for (size_t j = 0; j < qbts.size(); j++) {
                        tuple<int, int> qbt = qbts[j];
                        (*in_progress)(get<0>(qbt), get<1>(qbt)) = (*in_progress)(get<0>(qbt), get<1>(qbt)) + 1;
                    }
                }

                for (size_t j = 0; j < chains.size(); j++) {
                    vector<tuple<int, int>> chain = chains[j];
                    gen.route_qbt(chain[chain.size() - 1]);
                    add_chain(chain, gen.get_key());
                }
            }      
            // (*ancillas).print(std::cout);
            // cout << endl;   
        }

          
        perform_bell_measurement();
        perform_syndrome_measurements();

        for (size_t i = 0; i < generators.size(); i++) {
            Generator gen = generators[i];
            if (gen.is_done()) {
                tuple<int, int> dest = gen.get_dest();
                (*dests)(get<0>(dest), get<1>(dest)) = false;

                vector<tuple<int, int>> qbts = gen.get_qbts();
                for (size_t j = 0; j < qbts.size(); j++) {
                    tuple<int, int> qbt = qbts[j];
                    (*in_progress)(get<0>(qbt), get<1>(qbt)) = (*in_progress)(get<0>(qbt), get<1>(qbt)) - 1;
                }
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


int Grid::route_independent_sets(vector<vector<tuple<int, int>>> gens) {
    // split the generators up into indepenent sets to be routed separately.
    // each independent set can be routed in one round. Total time is thus number of independent sets

    size_t orig_size = gens.size();
    for (size_t rounds = 0; rounds < (orig_size + 1); rounds++) {
        vector<vector<tuple<int, int>>> ind_set = {gens[0]};
        gens.erase(begin(gens));

        for (auto it = begin(gens); it != end(gens);) {
            bool add = true;

            for (auto it2 = begin(ind_set); it2 != end(ind_set); it2++) {
                vector<tuple<int, int>> intersection;
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