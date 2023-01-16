#ifndef point_H
#define point_H

#include <tuple>
#include <cmath>

using namespace std;

class Point {
    public:
        int x, y;
        Point(int inp_x, int inp_y);
        Point();
        bool operator<(const Point& p) const;
        bool operator==(const Point& p) const;
        bool operator!=(const Point& p) const;      
        double distance(const Point& p) const;  
};

#endif