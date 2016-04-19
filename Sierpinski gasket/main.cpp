#include <iostream>
#include "bitmap_image.hpp"
#include <stdexcept>
#include <cmath>

using namespace std;

const bool DEBUG = false;

const double pi = atan(1)*4;
const double COS_1st = cos(-pi/2.0 + 2.0*pi/3.0);
const double SIN_1st = sin(-pi/2.0 + 2.0*pi/3.0);
const double COS_2nd = cos(-pi/2.0 + 4.0*pi/3.0);
const double SIN_2nd = sin(-pi/2.0 + 4.0*pi/3.0);
const double TAN_PI_BY_6 = tan(pi/6.0);
const double SIN_PI_BY_6 = sin(pi/6.0);

class triangle
{
public:
    const double x;
    const double y;
    const double radius;
    const int orientation;

    triangle () :x(0), y(0), radius(0), orientation(0) {};
    triangle (double x, double y, const double r, const int o);
    void draw(bitmap_image& image, int r, int g, int b);
};

triangle::triangle (double x, double y, const double r, const int o)
 :x(x), y(y), radius(r), orientation(o)
{
    if (orientation!=-1 && orientation!=1) throw runtime_error{"invalid orientation"};
}

void triangle::draw(bitmap_image& image, int r, int g, int b)
{
    if (DEBUG) {
        image_drawer draw(image);
        draw.pen_width(3);
        draw.pen_color(255,0,0);
        draw.circle(x,y,6);
    }

    // set all pixels within this triangle to given colour
    int64_t cy,cx; // current x,y coord
    int64_t start_y,end_y;
    switch (orientation) {
    case 1:
        start_y = y - radius;
        end_y = y + radius*SIN_PI_BY_6;
        break;
    case -1:
        start_y = y - radius*SIN_PI_BY_6;
        end_y = y + radius;
        break;
    default:
        throw runtime_error{"invalid orientation"};
    }

    for (cy=start_y; cy<end_y; ++cy) {
        int64_t hx; // half the length of the horizontal line connecting the vertical line through (x,y) with the sides of the triangle
        // calculate cx
        switch (orientation) {
        case 1:
            hx = (cy-(y-radius))*TAN_PI_BY_6;
            break;
        case -1:
            hx = (y+radius-cy)*TAN_PI_BY_6;
            break;
        default:
            throw runtime_error{"invalid orientation"};
        }

        // write the pixels
        for (cx=x-hx; cx<x+hx; ++cx)
            image.set_pixel(cx,cy, r,g,b);
    }
}

class gasket_t
{
public:
    vector<vector<triangle>> level; // level[n] = the 3^n triangles in level n of the gasket
    const double min_radius;
    gasket_t (const double origin_x, const double origin_y, const double radius, const double min_radius);
};

gasket_t::gasket_t (const double origin_x, const double origin_y, const double radius, const double min_radius)
 : min_radius(min_radius)
{
    // calculate levels of the gasket up until the radius of the given level is <1
    level.resize(1);
    level[0].emplace_back(origin_x, origin_y, radius, 1);

    double r = radius; // of current level
    while (true) {
        r = r/2.0;
        if (r<min_radius) break;

        const int64_t n = level.size() - 1;
        level.resize(level.size()+1);
        for(triangle& t : level[n]) {
            // add three new triangles to next level
            level[n+1].emplace_back(t.x, t.y-r, r, 1);
            level[n+1].emplace_back(t.x+r*COS_1st, t.y+r*SIN_1st, r, 1);
            level[n+1].emplace_back(t.x+r*COS_2nd, t.y+r*SIN_2nd, r, 1);
        }
    }

    cout << "constructed gasket with " << level.size() << " levels" << endl;
}

int main()
{
    // image dimensions
    const int64_t x_res = 2*10236; // pixels
    const int64_t y_res = 2*7087;

    // gasket dimensions
    const double x_origin = (double)x_res/2.0;
    const double y_origin = (double)y_res/1.6;
    const double radius = (double)min(x_res,y_res)*0.6;
    const double min_radius = 20;

    const double offset = 1.0; // artificial reduction of hole radius for visual effect

    // construct gasket, orientated upwards
    gasket_t gasket {x_origin, y_origin, radius, min_radius};
    bitmap_image image {x_res, y_res};
    image.set_all_channels(255,255,255);

    // draw level 0 in black
    gasket.level[0].front().draw(image, 0,0,0);

    // punch out the holes in white
    for(vector<triangle>& level : gasket.level) {
        for (triangle& t : level) {
            triangle hole {t.x, t.y, t.radius/2.0-offset, -1};
            hole.draw(image, 255,255,255);
        }
    }

    image.save_image("sierpinski_gasket.bmp");

}
