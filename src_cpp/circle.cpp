#include "circle.h"

Circle::Circle(vector<Point> points) {
    Circle c = Circle(0,0,-1);

    for (size_t i = 0; i < points.size(); i++) {
        if (c.r == -1 || !is_in_circle(c, points[i])) {
            c = make_circle_one_point(vector<Point>(points.begin(), points.begin()+i+1), points[i]);
        }
    }

    cx = c.cx;
    cy = c.cy;
    r = c.r;
}

Circle::Circle(double inp_cx, double inp_cy, double inp_r) :
    cx(inp_cx), cy(inp_cy), r(inp_r) {}

Circle::Circle() : cx(-1), cy(-1), r(-1) {}

bool Circle::operator== (const Circle& c) const {
    return make_tuple(cx, cy, r) == make_tuple(c.cx, c.cy, c.r);
}

bool Circle::operator!= (const Circle& c) const {
    return !Circle::operator==(c);
}

Circle Circle::make_circle_one_point(vector<Point> points, Point p) {
    Circle c = Circle(p.x, p.y, 0.0);

    for (size_t i = 0; i < points.size(); i++) {
        if (!is_in_circle(c, points[i])) {
            if (c.r == 0.0) {
                c = make_diameter(p, points[i]);
            } else {
                c = make_circle_two_points(vector<Point>(points.begin(), points.begin()+i+1), 
                    p, points[i]);
            }
        }
    }

    return c;
}

Circle Circle::make_circle_two_points(vector<Point> points, Point p, Point q) {
    Circle circ = make_diameter(p, q);
    Circle left = Circle(0,0,-1);
    Circle right = Circle(0,0,-1);
    int px = p.x, py = p.y;
    int qx = q.x, qy = q.y;

    for (size_t i = 0; i < points.size(); i++) {
        Point r = points[i];
        if (is_in_circle(circ, r)) continue;

        double cross = cross_product(px, py, qx, qy, r.x, r.y);
        Circle c = make_circumcircle(p, q, r);
        if (c.r == -1) {
            continue;
        } else if (cross > 0.0 && (left.r == -1 || (cross_product(px, py, qx, qy, c.cx, c.cy) > cross_product(px, py, qx, qy, left.cx, left.cy)))) {
            left = c;
        } else if (cross < 0.0 && (right.r == -1 || (cross_product(px, py, qx, qy, c.cx, c.cy) < cross_product(px, py, qx, qy, right.cx, right.cy)))) {
            right = c;
        }
    }

    if (left.r == -1 && right.r == -1) {
        return circ;
    } else if (left.r == -1) {
        return right;
    } else if (right.r == -1) {
        return left;
    } else {
        return (left.r <= right.r) ? left : right;
    }
}

Circle Circle::make_diameter(Point a, Point b) {
    double cx = (a.x + b.x) / 2;
    double cy = (a.y + b.y) / 2;

    double r0 = hypot(cx - a.x, cy - a.y);
    double r1 = hypot(cx - b.x, cy - b.y);
    
    return Circle(cx, cy, max(r0, r1));
}

Circle Circle::make_circumcircle(Point a, Point b, Point c) {
    double ox = (min({a.x, b.x, c.x}) + max({a.x, b.x, c.x})) / 2;
    double oy = (min({a.y, b.y, c.y}) + max({a.y, b.y, c.y})) / 2;
    
    double ax = a.x - ox;
    double ay = a.y - oy;
    double bx = b.x - ox;
    double by = b.y - oy;
    double cx = c.x - ox;
    double cy = c.y - oy;

    double d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2.0;
	if (!d) return Circle(0,0,-1);
    double x = ox + ((ax*ax + ay*ay) * (by - cy) + (bx*bx + by*by) * (cy - ay) + (cx*cx + cy*cy) * (ay - by)) / d;
	double y = oy + ((ax*ax + ay*ay) * (cx - bx) + (bx*bx + by*by) * (ax - cx) + (cx*cx + cy*cy) * (bx - ax)) / d;
    double ra = hypot(x - a.x, y - a.y);
    double rb = hypot(x - b.x, y - b.y);
    double rc = hypot(x - c.x, y - c.y);

    return Circle(x, y, max({ra, rb, rc}));
}

bool Circle::is_in_circle(Circle c, Point p) {
    // assuming c is not NULL
    return hypot(p.x - c.cx, p.y - c.cy) <= (c.r * (1 + 1e-14));
}

double Circle::cross_product(double x0, double y0, double x1, double y1, double x2, double y2) {
    return ((x1 - x0) * (y2 - y0)) - ((y1 - y0) * (x2 - x0));
}
