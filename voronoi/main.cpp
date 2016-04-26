#include <iostream>
#include <random>
#include <iomanip>
#include "bitmap_image.hpp"

using namespace std;

const int64_t x_res = 10236;
const int64_t y_res = 7087;

bitmap_image image {x_res, y_res};

const int64_t Epoints = 500;

default_random_engine random_engine {5}; // seed
uniform_real_distribution<double> std_uniform_dist {0.0,1.0};

int sample_Poisson (double mean)
{
    poisson_distribution<int> std_poisson_dist {mean};
    return std_poisson_dist(random_engine);
}

double sample_uniform ()
{
    return std_uniform_dist(random_engine);
}

class point_t
{
public:
    const int x;
    const int y;
    unsigned char colour;
    point_t (const int x, const int y) : x(x),y(y) {colour = 60 + int(sample_uniform()*20)*8;};
};

int main()
{
	cout << setprecision(2) << fixed;
    image.set_all_channels(255,255,255);

    vector<point_t> PPP;

    // construct the PPP
    const int n_point = sample_Poisson(Epoints);
    int i;
    for (i=0; i<n_point; ++i)
    {
        int x = sample_uniform() * x_res;
        int y = sample_uniform() * y_res;
        PPP.emplace_back(x, y);
    }
    cout << "PPP has " << n_point << " points" << endl;

    // colour the cells
    int x,y;
    for (x=0; x<x_res; ++x) {
    for (y=0; y<y_res; ++y) {
        int nearest = -1;
        double nearest_sqr_distance = 1e15;
        for (i=0; i<PPP.size(); ++i) {
            point_t& p = PPP[i];
            double d = (p.x-x)*(p.x-x) + (p.y-y)*(p.y-y);
            if (d<nearest_sqr_distance) {
                nearest = i;
                nearest_sqr_distance = d;
            }
        }
        unsigned char c = PPP[nearest].colour;
        image.set_pixel(x,y, c,c,c);
     }
    if (x%(x_res/100)==0) cout << "\rcoloured " << (double)x/(double)x_res;
        }
    cout << "\rcoloured " << 1.00;
    cout << endl;
    cout << "done" << endl;

    image.save_image("voronoi.bmp");



}
