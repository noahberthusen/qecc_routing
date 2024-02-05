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
        Grid(int M, int k);
        ~Grid();
        vector<Point> find_chain(const mat<int>& grid, Point site1, Point site2);
        void add_chain(vector<Point> chain, int gen);
        static void tmp_add_chain(mat<int>& grid, vector<Point> chain);
        void perform_bell_measurement();
        void perform_syndrome_measurements();
        vector<vector<Point>> route_generator(vector<Point> gen, Point prior_dest=Point());
        int greedy_route_set(vector<vector<Point>> gens);
        int route_independent_sets(vector<vector<Point>> gens);
        const mat<int>& get_ancillas();

        // functions for adding generators back in
        void greedy_route_reccuring(vector<vector<Point>> gens, int max_rounds, int restart);
        void reset_generator(Generator& gen);

    private:
        int M;
        int k;

        vector<Generator> generators;
        vector<Generator> finished_gens;
        vector<vector<vector<Point>>> full_chains;
        vector<vector<tuple<Point, Point>>> bell_pairs;

        mat<int>* ancillas;
        mat<bool>* dests;
        mat<int>* in_progress;

};


#endif