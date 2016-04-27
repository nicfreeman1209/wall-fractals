#include <iostream>
#include <deque>
#include "bitmap_image.hpp"
#include "random_sampler.hpp"

using namespace std;

const int64_t x_res = 2*10236; // pixels
const int64_t y_res = 2*7087;
bitmap_image image {x_res, y_res};

sampler_t sample {138};

const double p = 1.0; // off->on rate per pixel
const double f = 0.00002; // fire rate per pixel, SoC should occur with p>>f
const int64_t n_pixels = x_res*y_res;

const double global_fire_rate = n_pixels*f;
const double e = 2.7182818284;

const double max_fires = 1e5;

class point_t
{
public:
	const int x;
	const int y;
	point_t (const int x, const int y) : x(x), y(y) {};
};

deque<point_t> fire;

inline void handle_neighbour (const int x, const int y)
{
	// add to fire, if on
	static unsigned char r,g,b;
	if (x>=0 && x<x_res && y>=0 && y<y_res) {
		image.get_pixel(x,y, r,g,b);
		if (r==255) fire.emplace_back(x,y);
	}
}

int main()
{
	int i,x,y;
	unsigned char r,g,b;

	cout << "setting initial state...";
	image.set_all_channels(0,0,0);
	for (x=0; x<x_res; ++x) {
	for (y=0; y<y_res; ++y) {
		if (sample.unif_real_01()>0.5) image.set_pixel(x,y, 255,255,255);
	}
	}
	cout << "done" << endl;

    int64_t n_fire = 0;
    double model_time = 0.0;

	double burn_time;
	int64_t switched_on;

	while (true)
	{
		// get time to next fire
		burn_time = sample.exp_dist(global_fire_rate);

		// switch on clocks in between now and the next fire
		switched_on = n_pixels * (1-pow(e,-p*burn_time)); // expected number of switched on clocks, use as proxy for real (p>>f so SLLN says its a good estimate)
		for (i=0; i<switched_on; ++i) {
			x = sample.unif_real_01() * x_res;
			y = sample.unif_real_01() * y_res;
			image.set_pixel(x,y, 255,255,255);
		}

		// burn
		int sx = sample.unif_real_01() * x_res;
		int sy = sample.unif_real_01() * y_res;
		fire.emplace_back(sx, sy);
		while (fire.size()>0)
		{
			point_t& p = fire.front();
			image.get_pixel(p.x,p.y, r,g,b);
			if (r==255) {
				handle_neighbour(p.x+1,p.y);
				handle_neighbour(p.x-1,p.y);
				handle_neighbour(p.x,p.y+1);
				handle_neighbour(p.x,p.y-1);
				image.set_pixel(p.x,p.y, 0,0,0);
			}
			fire.pop_front();
		}

		// book-keeping
		model_time += burn_time;
		++n_fire;
		if (n_fire%1000==0) cout << "\rn_fire = " << n_fire;
		if (n_fire>=max_fires) break;
	}
	cout << "\rn_fire = " << n_fire << endl;
	cout << "model_time = " << model_time << endl;


	image.save_image("drossel-schwabl.bmp");
	return 0;
}
