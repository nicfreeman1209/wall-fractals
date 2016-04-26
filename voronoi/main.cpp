#include <iostream>
#include <random>
#include <iomanip>
#include <set>
#include <vector>
#include <stdexcept>
#include "bitmap_image.hpp"

using namespace std;

const int64_t x_res = 2*10236;
const int64_t y_res = 2*7087;

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
    const int id;
    int colour;
    point_t (const int x, const int y, const int id) : x(x),y(y),id(id) {colour = -1;};
};

int main()
{
	cout << setprecision(2) << fixed;
    image.set_all_channels(255,255,255);

    vector<point_t> PPP;
    vector<set<int>> graph; // bidirectional

    // construct the PPP
    const unsigned int n_point = sample_Poisson(Epoints);
    unsigned int i,j,k;
    for (i=0; i<n_point; ++i)
    {
        int x = sample_uniform() * x_res;
        int y = sample_uniform() * y_res;
        PPP.emplace_back(x, y, i);
    }
    cout << "PPP has " << PPP.size() << " points" << endl;
    graph.resize(PPP.size());

    // construct the graph of neighbours
    // two point of the PPP having neighbouring Voronoi cells iff there exists a point equidistant from both for which there is no other point closer
    // we approximate and hope
    const double dt_step = 1;
    const double dt_max = 500;

    int n_edges = 0;
    for (i=0; i<PPP.size(); ++i) {
    point_t& p1 = PPP[i];
    for (j=0; j<i; ++j) {
        point_t& p2 = PPP[j];
        // vector from p1 to p2
        double vx = p1.x - p2.x;
        double vy = p1.y - p2.y;
        // tangent normal
        double tx = -vy / sqrt(vx*vx+vy*vy);
        double ty = vx / sqrt(vx*vx+vy*vy);
        // midpoint
        double midx = (p1.x+p2.x)/2.0;
        double midy = (p1.y+p2.y)/2.0;

        // walk along the line defined midpoint + dt_cur * tangent_normal
        double dt_cur = -dt_max;
        bool neighbours;
        while (dt_cur < dt_max) {
            neighbours = true;

            double mx = midx + dt_cur*tx;
            double my = midy + dt_cur*ty;
            double d12_sqr = (mx-p1.x)*(mx-p1.x) + (my-p1.y)*(my-p1.y);
            for (k=0; k<PPP.size(); ++k) {
                if (k==i || k==j) continue;
                point_t& p3 = PPP[k];
                double d_sqr = (mx-p3.x)*(mx-p3.x) + (my-p3.y)*(my-p3.y);
                if (d_sqr<d12_sqr) {
                    neighbours = false;
                    break;
                }
            }
            if (neighbours) break;
            dt_cur += dt_step;
        }
        if (neighbours) {
            graph[i].insert(j);
            graph[j].insert(i);
            ++n_edges;
        }
    }
    }
    cout << "found " << n_edges << " edges" << endl;

    // colour the vertices
    // avoid giving neighbours the same colour
    int n_colour = 0;
    for (point_t& p : PPP) {
        set<int> neighbour_colours;
        for (int n_id : graph[p.id]) {
            point_t& n = PPP[n_id];
            neighbour_colours.insert(n.colour);
        }
        if (neighbour_colours.size()==0) throw runtime_error{"found isolated vertex"};
        // assign the lowest int that is not already taken by a neighbour
        int colour = 0;
        while (neighbour_colours.find(colour)!=neighbour_colours.end()) {
            ++colour;
        }
        p.colour = colour;
        n_colour = max(n_colour, colour);
    }
    ++n_colour;
    cout << "coloured PPP with " << n_colour << " colours" << endl;

    // construct random permutation of length n_colour
    vector<int> colour_map;
    colour_map.resize(n_colour);
    for (i=0; i<colour_map.size(); ++i) colour_map[i] = i;
    for (i=0; i<colour_map.size()-1; ++i) {
        int j = i + sample_uniform()*(colour_map.size()-i);
        int temp = colour_map[i];
        colour_map[i] = colour_map[j];
        colour_map[j] = temp;
    }

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
        unsigned char c = 50 + colour_map[PPP[nearest].colour] * (175.0/n_colour);
        image.set_pixel(x,y, c,c,c);
    }
    if (x%(x_res/100)==0) cout << "\rcoloured " << (double)x/(double)x_res;
    }
    cout << "\rcoloured " << 1.00;
    cout << endl;
    cout << "done" << endl;

    // draw edges (for debugging)
    /*
    image_drawer draw {image};
    draw.pen_color(255,0,0);
    draw.pen_width(3);
    for (i=0; i<graph.size(); ++i) {
        point_t& p1 = PPP[i];
        for (int j : graph[i]) {
            point_t& p2 = PPP[j];
            draw.line_segment(p1.x,p1.y, p2.x,p2.y);
        }
    }
    */

    image.save_image("voronoi.bmp");
}
