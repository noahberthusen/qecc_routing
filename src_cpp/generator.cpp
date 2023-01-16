#include "generator.h"

Generator::Generator(vector<Point> qbts, int i) {
    key = i;
    dest = Point();

    all_qbts = qbts;
    qbts_to_route = qbts;
    routed = {};
}

Generator::Generator(const Generator& gen) {
    key = gen.key;
    dest = gen.dest;

    all_qbts = gen.all_qbts;
    qbts_to_route = gen.qbts_to_route;
    routed = gen.routed;
}

Generator::~Generator() {}

void Generator::route_qbt(Point qbt) {
    vector<Point>::iterator pos = find(begin(qbts_to_route), end(qbts_to_route), qbt);
    if (pos != end(qbts_to_route)) qbts_to_route.erase(pos);
    routed.push_back(qbt);
}

int Generator::get_key() const {
    return key;
}

Point Generator::get_dest() const {
    return dest;
}

void Generator::set_dest(Point new_dest) {
    dest = new_dest;
}

bool Generator::is_done() const {
    return (qbts_to_route.size() == 0);
}

size_t Generator::num_routed() const {
    return routed.size();
}

vector<Point> Generator::get_qbts() const {
    return all_qbts;
}

vector<Point> Generator::get_routed_qbts() const {
    return routed;
}

vector<Point> Generator::get_qbts_to_route() const {
    return qbts_to_route;
}