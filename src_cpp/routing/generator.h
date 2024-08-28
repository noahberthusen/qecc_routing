#ifndef generator_H
#define generator_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include "point.h"
#include "circle.h"

using namespace std;

class Generator {
    public:
        int key, start, cycle;
        Point dest;
        vector<Point> all_qbts;
        vector<Point> qbts_to_route;
        vector<Point> routed;
        Circle c;

        Generator(vector<Point> qbts, int i);
        Generator(const Generator& gen);
        ~Generator();
        void route_qbt(Point qbt);
        size_t num_routed() const;
        bool is_done() const;
        void reset();
};

#endif