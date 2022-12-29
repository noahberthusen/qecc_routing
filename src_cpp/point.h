#ifndef point_H
#define point_H

#include <tuple>

using namespace std;

class Point {
    public:
        int x, y;
        Point(int inp_x, int inp_y);
        Point();
        bool operator< (const Point& pointObj) const;
        bool operator== (const Point& pointObj) const;
        bool operator!= (const Point& pointObj) const;        
};

#endif