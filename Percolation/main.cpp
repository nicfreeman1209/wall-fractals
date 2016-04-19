#include <iostream>
#include "bitmap_image.hpp"
#include "random_sampler.hpp"
#include <stdexcept>
#include <cmath>
#include <deque>
#include <map>

// image dimensions
const int64_t x_res = 2*10236; // pixels
const int64_t y_res = 2*7087;
bitmap_image image {x_res, y_res};
bitmap_image percolation {x_res, y_res};

const double p_percolation = 0.5924; //0.592746 is approx threshold, go slightly under to avoid large clusters all touching each other

class point_t {
public:
    const int64_t x;
    const int64_t y;
    point_t (const int64_t x, const int64_t y) :x(x),y(y) {};
};

bool is_within_image(const int x, const int y)
{
    return (0<=x && x<x_res && 0<=y && y<y_res);
}

vector<string> colours = {
    // there is no better method of colouring than a hard coded list
    "FF0000","00AA00","0000FF","000000","0076FF","010067","95003A","007DB5",
    "FF00F6","FFEEE8","774D00","90FB92","01FFFE","D5FF00","FF937E","6A826C",
    "FF029D","FE8900","7A4782","7E2DD2","85A900","FF0056","A42400","00AE7E",
    "683D3B","BDC6FF","263400","BDD393","00B917","9E008E","001544","C28C9F",
    "FF74A3","01D0FF","004754","E56FFE","788231","0E4CA1","91D0CB","BE9970",
    "968AE8","BB8800","43002C","DEFF74","00FFC6","FFE502","620E00","008F9C",
    "98FF52","7544B1","B500FF","00FF78","FF6E41","005F39","6B6882","5FAD4E",
    "A75740","A5FFD2","FFB167","009BFF","E85EBE"
};

vector<unsigned char> RGB_colour (const unsigned int colour_id)
{
    if (colour_id<0 || colour_id>=colours.size()) throw runtime_error{"invalid colour id"};

    int num = stoi(colours[colour_id], 0, 16); // hex
    int r = num / 0x10000; // 0..255
    int g = (num / 0x100) % 0x100;
    int b = num % 0x100;

    vector<unsigned char> colour;
    colour.resize(3);
    colour[0] = r;
    colour[1] = g;
    colour[2] = b;
    return colour;
}

bool colour_cluster (point_t& c, vector<unsigned char>& col)
{
    unsigned char r,g,b;
    percolation.get_pixel(c.x,c.y, r,g,b);
    if (r==0) return false; // closed site -> no cluster
    if (g==1) return true; // already in a cluster

    // map out new cluster
    deque<point_t> q;
    q.push_back(c);
    image.set_pixel(c.x,c.y, col[0], col[1], col[2]);
    percolation.set_pixel(c.x,c.y, r,1,0);

    int i,j;
    while (q.size()>0)
    {
        point_t& p = q.front();
        for (i=p.x-1; i<=p.x+1; i+=1) {
        for (j=p.y-1; j<=p.y+1; j+=1) {
            if (i!=p.x && j!=p.y) continue;
            if (i==p.x && j==p.y) continue;
            if (!is_within_image(i,j)) continue;
            percolation.get_pixel(i,j, r,g,b);
            if (r==1 && g==0) {
                // open site, add to this cluster
                image.set_pixel(i,j, col[0],col[1],col[2]);
                percolation.set_pixel(i,j, 1,1,0);
                q.emplace_back(i,j);
            }
        }
        }
        q.pop_front();
    }
    return true;
}

bool check_cluster_size (point_t& c, const int wanted_cluster_size)
{
    // find out if drawing a cluster starting from the point c would add a new cluster with more points that wanted_cluster_size
    // don't change the red/green channels whilst doing so
    unsigned char r,g,b;
    percolation.get_pixel(c.x,c.y, r,g,b);
    if (r==0) return false; // closed site -> no cluster
    if (g==1) return false; // already in a cluster -> no new cluster

    // map out new cluster, using blue channel to mark used points
    // afterwards, clean out blue channel
    deque<point_t> q;
    q.push_back(c);
    percolation.set_pixel(c.x,c.y, r,g,1);
    int64_t n = 1;

    int i,j;
    while (q.size()>0)
    {
        point_t& p = q.front();
        for (i=p.x-1; i<=p.x+1; i+=1) {
        for (j=p.y-1; j<=p.y+1; j+=1) {
            if (i!=p.x && j!=p.y) continue;
            if (i==p.x && j==p.y) continue;
            if (!is_within_image(i,j)) continue;
            percolation.get_pixel(i,j, r,g,b);
            if (r==1 && g==0 && b==0) {
                // open site, add to this cluster
                percolation.set_pixel(i,j, 1,g,1);
                ++n;
                if (n>=wanted_cluster_size) return true;
                q.emplace_back(i,j);
            }
        }
        }
        q.pop_front();
    }

    // remove marks
    int x,y;
    for (x=0; x<x_res; ++x) {
    for (y=0; y<y_res; ++y) {
        percolation.get_pixel(x,y, r,g,b);
        percolation.set_pixel(x,y, r,g,0);
    }
    }
    return false;
}

int main()
{
    sampler_t sample {138};

    // image
    image.set_all_channels(255,255,255);
    percolation.set_all_channels(0,0,0);

    // perform percolation, record into r channel as 0 or 1
    int x,y;
    cout << "calculating percolation...";
    for (x=0; x<x_res; ++x) {
    for (y=0; y<y_res; ++y) {
        unsigned char r = (sample.unif_real_01()<p_percolation) ? 1 : 0;
        percolation.set_pixel(x,y, r,0,0);
    }
    }
    cout << "done" << endl;

    // specially coloured clusters
    const int special_clusters = 4;

    // set a few coloured pixels as seeds
    // make sure the correspond to large-ish clusters near the centre of the image
    int i;
    const int max_attempts = 50;
    int attempts = 0;
    const int wanted_cluster_size = 1e5;
    cout << "colouring special clusters...";
    for (i=0; i<special_clusters; ++i) {
        while (true) {
            x = x_res/10 + 8*sample.unif_real_01()*x_res/10;
            y = y_res/10 + 8*sample.unif_real_01()*y_res/10;
            point_t c {x,y};
            if (check_cluster_size(c, wanted_cluster_size)) break;
            ++attempts;
            if (attempts>=max_attempts) {
                cout << "warning: could not find large enough cluster (i=" << i << ")" << endl;
                break;
            }

        }

        point_t p {x,y};
        vector<unsigned char> col = RGB_colour(i);

        percolation.set_pixel(p.x,p.y, 1,0,0); // force on!
        colour_cluster(p, col);
        cout << ".";
    }
    cout << "done" << endl;

    image.save_image("percolation.bmp");

    // colour the rest in shades of light grey
    cout << "colouring other clusters...";
    for (x=0; x<x_res; ++x) {
    for (y=0; y<y_res; ++y) {
        unsigned char grey = 150 + sample.unif_real_01() * 75;
        point_t p {x,y};
        vector<unsigned char> col = {grey, grey, grey};
        colour_cluster(p, col);
    }
    if (x%(x_res/20)==0) cout << ".";
    }
    cout << "done" << endl;

    image.save_image("percolation.bmp");
}
