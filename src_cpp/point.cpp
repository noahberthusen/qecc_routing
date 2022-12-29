#include "point.h"

Point::Point(int inp_x, int inp_y): x(inp_x), y(inp_y) {}

Point::Point(): x(-1), y(-1) {}

bool Point::operator< (const Point& pointObj) const {
    return make_tuple(y, x) < make_tuple(pointObj.y, pointObj.x);
}

bool Point::operator== (const Point& pointObj) const {
    return make_tuple(y, x) == make_tuple(pointObj.y, pointObj.x);
}

bool Point::operator!= (const Point& pointObj) const {
    return !Point::operator==(pointObj);
}