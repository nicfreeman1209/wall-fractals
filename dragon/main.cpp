#include <iostream>
#include "bitmap_image.hpp"
#include <cmath>
#include <list>

const int x_res = 2*10236;
const int y_res = 2*7087;

const int level = 16;

const double r = 24; // lattice size

// lazy by-hand positioning
const int o_x = 11*x_res/28;
const int o_y = 25*y_res/36;

bitmap_image image {x_res, y_res};
image_drawer draw {image};

using namespace std;

bool is_in_image(const int x, const int y)
{
    return (x>=0 && x<x_res && y>0 && y<y_res);
}

int decode_direction(int dir)
{
    if (dir%2==0) return 0;
    return (dir-2);
}

void draw_line(int x1, int y1, int x2, int y2)
{
    if (is_in_image(x1,y1) && is_in_image(x2,y2))
        draw.line_segment(x1,y1, x2,y2);
}

int main()
{
    image.set_all_channels(255,255,255);
    draw.pen_color(0,0,0);
    draw.pen_width(3);

    // construct the sequence of turns
    list<int> pattern;

    pattern.push_back(1);
    int i;
    for (i=0; i<level; ++i) {
        list<int> swapped_pattern;
        for (int dir : pattern)
                swapped_pattern.push_back(1-dir);

        pattern.push_back(1);
        swapped_pattern.reverse();
        pattern.splice(pattern.end(), swapped_pattern);
    }
    cout << "pattern length: " << pattern.size() << "=2*" << (pattern.size()-1)/2 << "+1" << endl;

    // construct the lines
    double x1,y1,x2,y2,x3,y3;
    double x_dir = 2;
    double y_dir = 1;
    int dir;

    list<int>::iterator it = pattern.begin();

    // first point
    x1 = o_x;
    y1 = o_y;

    // second point
    dir = *it;
    ++it;

    if (dir==1) {
        // turn right
        x_dir += 1;
        if (x_dir==4) x_dir = 0;
        y_dir += 1;
        if (y_dir==4) y_dir = 0;
    }
    else {
        // turn left
        x_dir -= 1;
        if (x_dir==-1) x_dir = 3;
        y_dir -= 1;
        if (y_dir==-1) y_dir = 3;
    }

    x2 = x1 + r * decode_direction(x_dir);
    y2 = y1 + r * decode_direction(y_dir);

    // main loop
    // we draw 'per corner', which means we need to track three points at once
    while(true)
    {
        // get third point from first and second
        dir = *it;
        ++it;

        if (dir==1) {
            // turn right
            x_dir += 1;
            if (x_dir==4) x_dir = 0;
            y_dir += 1;
            if (y_dir==4) y_dir = 0;
        }
        else {
            // turn left
            x_dir -= 1;
            if (x_dir==-1) x_dir = 3;
            y_dir -= 1;
            if (y_dir==-1) y_dir = 3;
        }

        double x3 = x2 + r * decode_direction(x_dir);
        double y3 = y2 + r * decode_direction(y_dir);

        // get various midpoints used in drawing the 'turning the corner' lines
        double x12 = (x1+x2)/2.0;
        double y12 = (y1+y2)/2.0;
        double x122 = (x1+3*x2)/4.0;
        double y122 = (y1+3*y2)/4.0;
        double x223 = (3*x2+x3)/4.0;
        double y223 = (3*y2+y3)/4.0;
        double x23 = (x2+x3)/2.0;
        double y23 = (y2+y3)/2.0;

        // lines
        draw_line(x12,y12, x122,y122);
        draw_line(x122,y122, x223,y223);
        draw_line(x223,y223, x23,y23);

        // did we finish?
        if (it==pattern.end()) break;

        // set new 1st & 2nd points
        x1 = x2;
        y1 = y2;
        x2 = x3;
        y2 = y3;
    }

    image.save_image("dragon.bmp");
}
