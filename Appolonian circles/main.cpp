#include <iostream>
#include "bitmap_image.hpp"
#include <cmath>
#include <complex>

using namespace std;

const double epsilon = 0.01;
const double pi = acos(-1.0);

class circle_t {
public:
    const double x;
    const double y;
    const double r;
    const bool external; // is one of its parents the outer circle?
    const int64_t id;

    const vector<int64_t> parents;

    circle_t (double x, double y, double r, int64_t id, bool e) : x(x), y(y), r(r), id(id), external(e) {};
    circle_t (double x, double y, double r, int64_t id, bool e, vector<int64_t> p) : x(x), y(y), r(r), id(id), external(e), parents(p) {};

    bool is_tangent_too (const circle_t& c);
};

bool circle_t::is_tangent_too (const circle_t& c)
{
    // quick and dirty
    double dist_sqr = (x-c.x)*(x-c.x) + (y-c.y)*(y-c.y);
    double rad_dist_sqr = (r+c.r)*(r+c.r);
    return (abs(dist_sqr-rad_dist_sqr) < epsilon);
}

class gasket_t {
public:
    const double o_x;
    const double o_y;
    const double radius;
    const double min_radius; //smallest circle we record

    vector<circle_t> circle;

    gasket_t (double x, double y, double r, double min_r);
private:
    void emplace_back_internal_circle (circle_t& c1, circle_t& c2, circle_t& c3);
    void emplace_back_external_circle (circle_t& c1, circle_t& c2);

};

void gasket_t::emplace_back_internal_circle (circle_t& c1, circle_t& c2, circle_t& c3)
{
    // return the circle which is internally tangent to c1,c2,c3
    // assume c1,c2,c3 are all externally tangent
    double k1 = 1.0 / c1.r;
    double k2 = 1.0 / c2.r;
    double k3 = 1.0 / c3.r;
    complex<double> C1 {c1.x, c1.y};
    complex<double> C2 {c2.x, c2.y};
    complex<double> C3 {c3.x, c3.y};
    double k4 = k1 + k2 + k3 + 2.0*sqrt(k1*k2+k2*k3+k1*k3);
    complex<double> C4 = (C1*k1 + C2*k2 + C3*k3 + 2.0*sqrt(k1*k2*C1*C2 + k2*k3*C2*C3 + k1*k3*C1*C3)) / k4;

    vector<int64_t> parents = {c1.id, c2.id, c3.id};
    circle.emplace_back(C4.real(), C4.imag(), 1.0/k4, circle.size(), false, parents);
}

void gasket_t::emplace_back_external_circle (circle_t& c1, circle_t& c2)
{
    // return the circle which is internally tangent to the outer circle and externally tangent to c1 and c2
    // assume c1,c2 are internally tangent to the outer circle and externally tangent to each other
    double k1 = -1.0 / circle[0].r;
    double k2 = 1.0 / c1.r;
    double k3 = 1.0 / c2.r;
    complex<double> C1 {circle[0].x, circle[0].y};
    complex<double> C2 {c1.x, c1.y};
    complex<double> C3 {c2.x, c2.y};
    double k4 = k1 + k2 + k3 + 2.0*sqrt(k1*k2+k2*k3+k1*k3);
    complex<double> C4 = (C1*k1 + C2*k2 + C3*k3 + 2.0*sqrt(k1*k2*C1*C2 + k2*k3*C2*C3 + k1*k3*C1*C3)) / k4;

    vector<int64_t> parents = {c1.id, c2.id};
    circle.emplace_back(C4.real(), C4.imag(), 1.0/k4, circle.size(), true, parents);
}

gasket_t::gasket_t (double o_x, double o_y, double radius, double min_radius)
 : o_x(o_x), o_y(o_y), radius(radius), min_radius(min_radius)
 {
    cout << "Apollonian gasket with center " << o_x << "," << o_y << " and radius " << radius << endl;

    // add the outermost circle
    circle.emplace_back(o_x,o_y,radius,0,false);

    // add the three initial internal circles (which have no parents)
    int64_t i;
        double rat = 2.0*sqrt(3.0)-3.0;
    for (i=0; i<3; ++i) {
        double r = radius*rat;
        double x = o_x + (radius-r) * cos(pi+2.0*pi*i/3.0); // don't remove the extra pi or something horrible goes wrong!
        double y = o_y + (radius-r) * sin(pi+2.0*pi*i/3.0);
        double id = i+1;
        circle.emplace_back(x,y,r,id,true);

    }

    // add the next four circles, which do have parents
    emplace_back_internal_circle(circle[1], circle[2], circle[3]);
    emplace_back_external_circle(circle[1], circle[2]);
    emplace_back_external_circle(circle[2], circle[3]);
    emplace_back_external_circle(circle[1], circle[3]);

    // iterate, starting from circle[4]
    // for each circle, add its three children on the end
    // stop when we reach circles smaller than min_radius
    int64_t n = 3;
    while (n<circle.size())
    {
        ++n;
        circle_t c = circle[n]; // have to take copies, because emplace invalidates iterators
        if (c.r<min_radius) continue;

        if (c.external)
        {
            if (c.parents.size()!=2) throw runtime_error{"invalid parents"};
            circle_t p0 = circle[c.parents[0]];
            circle_t p1 = circle[c.parents[1]];

            emplace_back_external_circle(c, p0);
            emplace_back_external_circle(c, p1);
            emplace_back_internal_circle(c, p0, p1);
        }
        else {
            if (c.parents.size()!=3) throw runtime_error{"invalid parents"};
            circle_t p0 = circle[c.parents[0]];
            circle_t p1 = circle[c.parents[1]];
            circle_t p2 = circle[c.parents[2]];

            emplace_back_internal_circle(c, p0, p1);
            emplace_back_internal_circle(c, p1, p2);
            emplace_back_internal_circle(c, p0, p2);
        }
    }

    cout << "Found " << circle.size() << " circles with min radius " << min_radius << endl;
 }

int main()
{
    const int64_t x_res = 2*10236;
    const int64_t y_res = 2*7087;

    const double x = (double)x_res/2.0;
    const double y = (double)y_res/2.0;
    const double r = (double)min(x_res,y_res)*0.45;

    const double min_radius = 1;

    gasket_t gasket {x,y,r , min_radius};

    bitmap_image image {x_res, y_res};
    image.set_all_channels(255,255,255);

    image_drawer draw {image};
    draw.pen_width(3);
    draw.pen_color(0,0,0);
    for (circle_t& c : gasket.circle) {
        //cout << "drawing: " << c.x << ", " << c.y << ", " << c.r << endl;
        draw.circle(c.x,c.y,c.r);
    }

    image.save_image("apollonian_gasket.bmp");


}
