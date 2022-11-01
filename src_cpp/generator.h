#ifndef generator_H
#define generator_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>

using namespace std;

class Generator {
    public:
        Generator(vector<tuple<int, int>> qbts, int i);
        Generator(const Generator& gen);
        ~Generator();
        vector<tuple<int, int>> get_qbts_to_route() const;
        void route_qbt(tuple<int, int> qbt);
        int get_key() const;
        void set_dest(tuple<int, int> new_dest);
        tuple<int, int> get_dest() const;
        size_t num_routed() const;
        bool is_done() const;
        vector<tuple<int, int>> get_qbts() const;
        vector<tuple<int, int>> get_routed_qbts() const;

    private:
        int key;
        tuple<int, int> dest;
        vector<tuple<int, int>> all_qbts;
        vector<tuple<int, int>> qbts_to_route;
        vector<tuple<int, int>> routed;
};

#endif