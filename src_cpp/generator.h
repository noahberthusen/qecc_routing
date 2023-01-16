#ifndef generator_H
#define generator_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include "point.h"

using namespace std;

class Generator {
    public:
        Generator(vector<Point> qbts, int i);
        Generator(const Generator& gen);
        ~Generator();
        vector<Point> get_qbts_to_route() const;
        void route_qbt(Point qbt);
        int get_key() const;
        void set_dest(Point new_dest);
        Point get_dest() const;
        size_t num_routed() const;
        bool is_done() const;
        vector<Point> get_qbts() const;
        vector<Point> get_routed_qbts() const;

    private:
        int key;
        Point dest;
        vector<Point> all_qbts;
        vector<Point> qbts_to_route;
        vector<Point> routed;
};

#endif