#ifndef grid_H
#define grid_H

#include <vector>
#include <tuple>
#include <deque>
#include <algorithm>
#include "mat.h"
#include "generator.h"

using namespace std;

class Grid {
    public:
        Grid(int N, int k);
        ~Grid();
        vector<tuple<int, int>> find_chain(const mat<int>& grid, tuple<int, int> site1, tuple<int, int> site2);
        void add_chain(vector<tuple<int, int>> chain, int gen);
        static void tmp_add_chain(mat<int>& grid, vector<tuple<int, int>> chain);
        void perform_bell_measurement();
        void perform_syndrome_measurements();
        vector<vector<tuple<int, int>>> route_generator(vector<tuple<int, int>> gen, tuple<int, int> prior_dest=make_tuple(-1, -1));
        int greedy_route_set(vector<vector<tuple<int, int>>> gens);
        int route_independent_sets(vector<vector<tuple<int, int>>> gens);
        const mat<int>& get_ancillas();

    private:
        int N;
        int k;

        vector<Generator> generators;
        vector<vector<vector<tuple<int, int>>>> full_chains;
        vector<vector<tuple<tuple<int, int>, tuple<int, int>>>> bell_pairs;

        mat<int>* ancillas;
        mat<bool>* dests;
        mat<int>* in_progress;

};


#endif