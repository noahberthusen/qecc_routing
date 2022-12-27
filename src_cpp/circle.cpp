#include "circle.h"

Circle::Circle(vector<tuple<int, int>> points) {

}

Circle::Circle(double cx, double cy, double r) {

}

double Circle::x() {
    
}


Circle Circle::make_diameter(tuple<int, int> a, tuple<int, int> b) {
    double cx = (get<0>(a) + get<0>(b)) / 2;
    double cy = (get<1>(a) + get<1>(b)) / 2;

    double r0 = hypot(cx - get<0>(a), cy - get<1>(a));
    double r1 = hypot(cx - get<0>(b), cy - get<1>(b));
    
    return Circle(cx, cy, max(r0, r1));

}

Circle Circle::make_circumcircle(tuple<int, int> a, tuple<int, int> b, tuple<int, int> c) {
    double ox = (min(get<0>(a), get<0>(b), get<0>(c)) + max(get<0>(a), get<0>(b), get<0>(c))) / 2;
    double oy = (min(get<1>(a), get<1>(b), get<1>(c)) + max(get<1>(a), get<1>(b), get<1>(c))) / 2;
    
    double ax = get<0>(a) - ox;
    double ay = get<1>(a) - oy;
    double bx = get<0>(b) - ox;
    double by = get<1>(b) - oy;
    double cx = get<0>(c) - ox;
    double cy = get<1>(c) - oy;

    double d = (ax * (by - cy)) * 2.0


}

bool Circle::is_in_circle(Circle c, tuple<int, int> p) {
    // assume c is not NULL
    return 0.0;
}

double Circle::cross_product(double x0, double y0, double x1, double y1, double x2, double y2) {
    return ((x1 - x0) * (y2 - y0)) - ((y1 - y0) * (x2 - x0));
}
