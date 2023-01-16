#include "point.h"

Point::Point(int inp_x, int inp_y): x(inp_x), y(inp_y) {}

Point::Point(): x(-1), y(-1) {}

bool Point::operator< (const Point& p) const {
    return make_tuple(x, y) < make_tuple(p.x, p.y);
}

bool Point::operator== (const Point& p) const {
    return make_tuple(x, y) == make_tuple(p.x, p.y);
}

bool Point::operator!= (const Point& p) const {
    return !Point::operator==(p);
}

double Point::distance(const Point& p) const {
    return sqrt(pow(x - p.x, 2) + pow(y - p.y, 2));
}