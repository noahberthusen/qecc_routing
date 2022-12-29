#ifndef circle_H
#define circle_H

#include "point.h"
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>

using namespace std;

class Circle {
    public:
        Circle(vector<Point> points);
        Circle(double cx, double cy, double r);
        static bool is_in_circle(Circle c, Point p);
        double cx, cy, r;

    private:
        Circle make_circle(vector<Point> points);
        Circle make_circle_one_point(vector<Point> points, Point p);
        Circle make_circle_two_points(vector<Point> points, Point p, Point q);
        Circle make_diameter(Point a, Point b);
        Circle make_circumcircle(Point a, Point b, Point c);
        double cross_product(double x0, double y0, double x1, double y1, double x2, double y2);
};


#endif