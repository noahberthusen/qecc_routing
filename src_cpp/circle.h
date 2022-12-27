#ifndef circle_H
#define circle_H

#include <vector>
#include <tuple>
#include <algorithm>

using namespace std;

class Circle {
    public:
        Circle(vector<tuple<int, int>> points);
        Circle(double cx, double cy, double r);
        
        double x();
        double y();
        double r();

    private:
        double x, y, r;
        Circle make_circle_one_point();

        Circle make_diameter(tuple<int, int> a, tuple<int, int> b);
        Circle make_circumcircle(tuple<int, int> a, tuple<int, int> b, tuple<int, int> c);
        bool is_in_circle(Circle c, tuple<int, int> p);
        double cross_product(double x0, double y0, double x1, double y1, double x2, double y2);
};


#endif