#include "generator.h"

Generator::Generator(vector<Point> qbts, int i) {
    key = i;
    dest = Point();
    cycle = 0;
    start = 0;

    all_qbts = qbts;
    qbts_to_route = qbts;
    routed = {};
    c = Circle(qbts);
}

Generator::Generator(const Generator& gen) {
    key = gen.key;
    dest = gen.dest;
    cycle = gen.cycle;

    all_qbts = gen.all_qbts;
    qbts_to_route = gen.qbts_to_route;
    routed = gen.routed;
    start = gen.start;
    c = gen.c;
}

Generator::~Generator() {}

void Generator::route_qbt(Point qbt) {
    vector<Point>::iterator pos = find(begin(qbts_to_route), end(qbts_to_route), qbt);
    if (pos != end(qbts_to_route)) qbts_to_route.erase(pos);
    routed.push_back(qbt);
}

bool Generator::is_done() const {
    return (qbts_to_route.size() == 0);
}

size_t Generator::num_routed() const {
    return routed.size();
}

void Generator::reset() {
    cycle++;
    dest = Point();

    qbts_to_route = all_qbts;
    routed = {};
}