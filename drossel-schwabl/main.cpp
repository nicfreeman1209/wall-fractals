#include <iostream>
#include <deque>
#include "bitmap_image.hpp"
#include "random_sampler.hpp"

using namespace std;

const bool write_video_frames = true;

const int64_t x_out = 300; // output resolution
const int64_t y_out = 200;
bitmap_image image_out {x_out, y_out};

const int64_t x_res = 12000; // internal resolution
const int64_t y_res = 8000;
bitmap_image image {x_res, y_res};

sampler_t sample {138};

const double p = 1.0; // off->on rate per pixel
const double f = 0.00001; // fire rate per pixel, SoC should occur with p>>f
const int64_t n_pixels = x_res*y_res;

const double global_fire_rate = n_pixels*f;
const double e = 2.7182818284;

const double max_fires = 5e4;

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
		if (r==0) fire.emplace_back(x,y);
	}
}

void write_image(bitmap_image& image_in, const string& file_name)
{
	// average pixels to requested size
	int64_t x,y;
	int64_t x_scale = x_res / x_out;
	int64_t y_scale = y_res / y_out;
	for (x=0; x<x_out; ++x)	{
	for (y=0; y<y_out; ++y) {
		int64_t val = 0.0;
		int64_t i,j;
		int64_t n=0;
		for (i=x*x_scale; i<(x+1)*x_scale; ++i) {
		for (j=y*y_scale; j<(y+1)*y_scale; ++j) {
			if (x>=x_res || y>=y_res) continue;
			unsigned char r,g,b;
			image_in.get_pixel(i,j, r,g,b);
			val += r;
			++n;
		}
		}
		val /= (double)n;
		val *= 255;
		image_out.set_pixel(x,y, val, val, val);
	}
	}

	image_out.save_image(file_name);
}

int main()
{
	int i,x,y;
	unsigned char r,g,b;

	cout << "setting initial state...";
	image.set_all_channels(255,255,255);
	for (x=0; x<x_res; ++x) {
	for (y=0; y<y_res; ++y) {
		if (sample.unif_real_01()>0.5) image.set_pixel(x,y, 0,0,0);
	}
	}
	cout << "done" << endl;

    int64_t n_fire = 0;
    double model_time = 0.0;

	double fps = 90.0;
	double frame_interval = 1.0/fps;
	double next_frame = 0.0;
	int64_t frame_number = 0;

	double time = 0.0;
	double burn_time;
	int64_t switched_on;

	cout << "\rcurrent frame: " << frame_number << " current time: " << time;
	while (true)
	{
		// record frame if needed
		while (next_frame <= time && write_video_frames)
		{
			write_image(image, "DS_" + to_string(frame_number) + ".bmp");
			++frame_number;
			next_frame += frame_interval;
		}
		cout << "\rcurrent frame: " << frame_number << " current time: " << time;

		// get time to next fire
		burn_time = sample.exp_dist(global_fire_rate);
		time += burn_time;

		// switch on clocks in between now and the next fire
		switched_on = n_pixels * (1-pow(e,-p*burn_time)); // expected number of switched on clocks, use as proxy for real (p>>f so SLLN says its a good estimate)
		for (i=0; i<switched_on; ++i) {
			x = sample.unif_real_01() * x_res;
			y = sample.unif_real_01() * y_res;
			image.set_pixel(x,y, 0,0,0);
		}

		// burn
		int sx = sample.unif_real_01() * x_res;
		int sy = sample.unif_real_01() * y_res;
		fire.emplace_back(sx, sy);
		while (fire.size()>0)
		{
			point_t& p = fire.front();
			image.get_pixel(p.x,p.y, r,g,b);
			if (r==0) {
				handle_neighbour(p.x+1,p.y);
				handle_neighbour(p.x-1,p.y);
				handle_neighbour(p.x,p.y+1);
				handle_neighbour(p.x,p.y-1);
				image.set_pixel(p.x,p.y, 255,255,255);
			}
			fire.pop_front();
		}

		// book-keeping
		model_time += burn_time;
		++n_fire;
		if (n_fire%1000==0) cout << "\rn_fire = " << n_fire;
		if (n_fire>=max_fires) break;

		if (time>=60.0) break;
	}
	cout << "\rn_fire = " << n_fire << endl;
	cout << "model_time = " << model_time << endl;


	image.save_image("drossel-schwabl.bmp");
	return 0;
}
