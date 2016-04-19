#include <iostream>
#include "bitmap_image.hpp"
#include "random_sampler.hpp"
#include <stdexcept>

using namespace std;

class BM_path {
public:
    const int64_t x_grid;
    const int64_t y_grid;
    const int64_t x_origin;
    const int64_t y_origin;
    const int64_t circle_radius;
    const int64_t circle_radius_squared;

    sampler_t sample;
    void walk ();
    bitmap_image* image;

    BM_path (const int64_t seed, bitmap_image* image);
private:
    void add_to_range(const int64_t, const int64_t);
};

BM_path::BM_path (const int64_t seed, bitmap_image* image)
 : x_grid(image->width()),
   y_grid(image->height()),
   x_origin(x_grid/2),
   y_origin(y_grid/2),
   circle_radius(0.45*min(x_grid,y_grid)),
   circle_radius_squared(circle_radius*circle_radius),
   image(image)
{
    sampler_t s {seed};
    sample = s;

    cout << "grid radii: " << x_grid << "," << y_grid << " with circle radius " << circle_radius << endl;
}

void BM_path::walk () {
    // start a BM at (0,0) and continue until it hits the disc of given radius
    int64_t x = x_origin;
    int64_t y = y_origin;
    int64_t iter = 0;

    while ((x-x_origin)*(x-x_origin) + (y-y_origin)*(y-y_origin) <= circle_radius_squared) {
        add_to_range(x,y);
        ++iter;
        if (iter%1000000==0) cout << "\r" << "iteration " << iter << ", displacement " << sqrt((x-x_origin)*(x-x_origin) + (y-y_origin)*(y-y_origin)) << "          ";

        double p = sample.unif_real_01();
        if (p<0.25)
            ++x;
        else if (p<0.5)
            --x;
        else if (p<0.75)
            --y;
        else
            ++y;
    }

    cout << "\r" << "hit circle at " << x-x_origin << "," << y-y_origin << "                          " << endl;
}

void BM_path::add_to_range(const int64_t x, const int64_t y)
{
    image->set_pixel(x,y,0,0,0);
}

int main()
{
    const int64_t x_res = 2*10236; // pixels
    const int64_t y_res = 2*7087;

    const int64_t seed = 516;

    bitmap_image image {x_res, y_res};
    image.set_all_channels(255,255,255);

    // draw the walk
    BM_path BM {seed, &image};
    BM.walk();

    // draw the circle
    image_drawer draw(image);
    draw.pen_width(3);
    draw.pen_color(0,0,0);
    draw.circle(BM.x_origin,BM.y_origin,BM.circle_radius);


    image.save_image("BM_to_disc.bmp");

    return 0;
}
