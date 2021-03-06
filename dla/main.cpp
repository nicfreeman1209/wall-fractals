#include <iostream>
#include <fstream>
#include <iomanip>
#include "bitmap_image.hpp"
#include <cmath>
#include <stdexcept>
#include <future>
#include <chrono>
#include <vector>
#include <list>
#include <random>
#include <cassert>

using namespace std;

// runs
const int64_t runs = 10;
const bool write_image_files = (runs==1) || false;

// image
const int64_t x_res = 1000*100; // pixels
const int64_t y_res = x_res;

// write output into tiles_sqrt^2 of equally sized tiles
// best to keep the tiles smaller than 30000x30000!
const int64_t tiles_sqrt = 4;
const int64_t tile_x_res = x_res / tiles_sqrt;
const int64_t tile_y_res = y_res / tiles_sqrt;

bitmap_image image {x_res, y_res}; // main data
bitmap_image tile {tile_x_res, tile_y_res}; // used only for writing tiles to file

// particles
const double particle_radius = 6; // values below 6 cause serious lattice effects; 0 -> particle==pixel

// jump distances
// for when a particle has enough space to make big circular jump without hitting the dla
// note: too big and book-keeping the regions costs too much, too small and the BMs cost too much + lattice effects occur
const double med_jump_radius = particle_radius + 5;
const double med_jump_radius_sqr = med_jump_radius*med_jump_radius;
const double long_jump_radius = 5*med_jump_radius;
const double long_jump_radius_sqr = long_jump_radius*long_jump_radius;

// border of image, in which we don't draw dla
const int64_t border = 9*particle_radius;

// dla/world radii
// these vary as the dla grows
double dla_radius;
double world_radius_sqr;
const double padding_factor = 10; // if a particles goes this far away (in multiples of current radii) from the centre of the circle inscribing dla, delete the particle (a "rejection")

// origin
const int64_t o_x = x_res/2;
const int64_t o_y = y_res/2;
const double o_x_dbl = o_x;
const double o_y_dbl = o_y;

// threading
// the simulation proceeds in stages:
//   during a stage each worker thread thread does its own calcs, with no knowledge of what the other threads do
//   at the end of each stage the main thread receives data back from all the worker threads, collates and processes it, then triggers the start of the next stage
// in more detail:
//   (1) during a stage, each thread simulates BMs until they hit the dla, and records the positions at which they hit
//       throughout the stage, the dla is *kept constant*, no new particles are added, and each worker thread uses the same picture of the dla
//   (2) at the end of the stage, each worker thread sends back a list of positions that were hit by the BMs
//       the main thread then goes through these positions, one by one, and adds corresponding particles to the dla (a "collate")
//       when the main thread finds a particle (call it X) trying to attach into a region of space that some other particle already been attached too (a "collision") it ignores X and moves on
//       after all positions are processed, the main thread initiates the next stage, and the worker threads now see the updated dla
//   (3) the reason for all this is that the dla (i.e. the bitmap_image, which is a C style array) can be read asynchronously, but cannot be written too asynchronously
//       so, this algorithm runs best on a single processor machine with lots of cores and lots of RAM
//   (4) see the real-time output for details of how many collisions/rejects are occurring
//
const int max_worker_threads = 10; // leave one thread free to bookkeep the worker threads (e.g. set this your number of logical cores, minus one)
const int max_BMs_per_thread_per_collate = 500; // hard limit
int BMs_per_thread_per_collate = -1;

// in progress info
const bool write_preview_images = false;
const int64_t preview_iter = 80*sqrt(x_res*y_res)/(1+particle_radius*particle_radius); // particles per preview (magic number)
const int64_t min_BMs_per_status_report = 1000;

// math consts
const double PI = atan(1)*4;
const int64_t random_seed = 1;


class dla_particle_t {
public:
	// a particle that is part of the dla (snapped to grid)
	const int x;
	const int y;
	dla_particle_t (const int x, const int y) :x(x), y(y) {};
};

class free_particle_t {
public:
	// a particle that is not (yet) part of the dla (not snapped to grid)
	double x;
	double y;
	bool hit_dla;
	free_particle_t (const double x, const double y) :x(x), y(y), hit_dla(false) {};

	double radial_displacement() {return sqrt((x-o_x)*(x-o_x)+(y-o_y)*(y-o_y));};
	void print() {cout << x << " " << y <<  " " << hit_dla << endl;};
	int64_t px () {return (int64_t)(x+0.5);};
	int64_t py () {return (int64_t)(y+0.5);};
};

inline bool is_in_world(free_particle_t& p, const int64_t r_sqr)
{
	return ((p.x-o_x)*(p.x-o_x) + (p.y-o_y)*(p.y-o_y) < r_sqr);
}

inline bool is_in_image(const int64_t x, const int64_t y)
{
	return (0<=x && x<x_res && 0<=y && y<y_res);
}

inline bool is_inside_dla_region(const int64_t x, const int64_t y)
{
	return (x>=border && x<x_res-border && y>=border && y<=y_res-border);
}

inline bool has_hit_dla(free_particle_t& p)
{
	// return true if (after being snapped to grid) we are within particle_radius of the boundary of an existing dla particle
	int64_t x = p.px();
	int64_t y = p.py();

	if (is_in_image(x,y)) {
		unsigned char r,g,b;
		image.get_pixel(x,y, r,g,b);
		return (g==127);
	}
	return false;
}

inline double unif_real_01 () {
	// thread safe rng!
	// constructs one rng per thread, has non-reproducible results because threading is non-deterministic
	uniform_real_distribution<double> distribution {0.0, 1.0};
	static mt19937 random_engine {random_seed};

	return distribution(random_engine);
}

free_particle_t BM_until_hit(double sample_dla_radius, double sample_world_radius_sqr)
{
	// start particle uniformly distributed on circle inscribing current dla
	const double sample_radius = sample_dla_radius + particle_radius + 1;
	const double sample_theta = unif_real_01()*2.0*PI;

	const double x = (double)o_x + sample_radius * cos(sample_theta);
	const double y = (double)o_y + sample_radius * sin(sample_theta);
	free_particle_t p {x,y};

	// key point: the exit distribution of BM started from the centre of a circle, is uniform
	double jump_radius;
	double theta;

	// now walk
	// return p with p.hit_dla=false if we exit world space
	// return p with p.hit_dla=true if p becomes neighbour to (current) dla
	while (true)
	{
		// give up if we exited the world
		if (!is_in_world(p, sample_world_radius_sqr)) return p;

		// stop if we hit the dla
		if (has_hit_dla(p)) {p.hit_dla=true; return p;};

		// if we are outside of the disc inscribing dla + long jump radius, do a 'kangaroo' circular jump
		// choose jump radius as large as possible
		if (p.radial_displacement()>dla_radius+particle_radius+long_jump_radius) {
			jump_radius = p.radial_displacement()-dla_radius;
			theta = unif_real_01()*2.0*PI;
			p.x += jump_radius * cos(theta);
			p.y += jump_radius * sin(theta);
			continue;
		}

		// if we are further than the long/med jump radius from dla, make a long/med jump
		unsigned char r = 0;
		unsigned char g = 0;
		unsigned char b = 0;

		int64_t px = p.px();
		int64_t py = p.py();
		if (is_in_image(px,py))
			image.get_pixel(px,py, r,g,b);

		if (b==0) {
			// long jump
			jump_radius = long_jump_radius;
			theta = unif_real_01()*2.0*PI;
			p.x += jump_radius * cos(theta);
			p.y += jump_radius * sin(theta);
			continue;
		}

		if (b==127) {
			// medium jump
			jump_radius = med_jump_radius;
			theta = unif_real_01()*2.0*PI;
			p.x += jump_radius * cos(theta);
			p.y += jump_radius * sin(theta);
			continue;
		}

		// advance BM by a single 'step'
		jump_radius = 1;
		theta = unif_real_01()*2.0*PI;
		p.x += jump_radius * cos(theta);
		p.y += jump_radius * sin(theta);
	}

	// should never happen
	cout << "error, particle stopped moving!" << endl;
	return p;
}

list<free_particle_t> thread_manager (double sample_dla_radius, double sample_world_radius_sqr)
{
	// run multiple BMs on the same thread (to mitigate the cost of async launches)
	int64_t i;
	list<free_particle_t> local_datapoint64_ts;
	for (i=0; i<BMs_per_thread_per_collate; ++i) {
		free_particle_t p =  BM_until_hit(sample_dla_radius, sample_world_radius_sqr);
		local_datapoint64_ts.push_back(p);
	}
	return local_datapoint64_ts;
}

void add_particle_to_dla(const int64_t x, const int64_t y)
{
	// use the green channel to mark the region in which a particle would (after this addition) be said to have hit the dla
	// use the blue channel to mark the region in which it is safe for particles to make med/long jumps
	// draw the particle in the red channel (so as the dla appears white in preview images)
	int64_t i,j;
	unsigned char r,g,b,new_g,new_b,new_r;
	for (i=x-2*particle_radius-long_jump_radius; i<=x+2*particle_radius+long_jump_radius; ++i) {
	for (j=y-2*particle_radius-long_jump_radius; j<=y+2*particle_radius+long_jump_radius; ++j) {
		if (!is_in_image(i,j)) continue;
		double dist = sqrt((x-i)*(x-i)+(y-j)*(y-j));
		if (dist>2*(double)particle_radius+long_jump_radius) continue;

		image.get_pixel(i,j, r,g,b);
		new_r = 0;
		new_b = 0;
		new_g = 0;

		// new rgbs reflecting only proximity to (x,y)
		// rgb values are monotonically increasing in proximity to the dla
		if (dist<=particle_radius) new_r=255; // -> inside the dla

		if (dist<=2*particle_radius+1) new_g=127; // -> a particle here would touch the particle we just added
		if (dist<=2*particle_radius) new_g=255; // -> a particle here would overlap the dla

		if (dist<=2*particle_radius+long_jump_radius) new_b=127; // -> we can't make a long jump here
		if (dist<=2*particle_radius+med_jump_radius) new_b=255; // -> we can't make a med jump here

		// merge with already known info
		if (r!=new_r || g!=new_g || b!=new_b)
			image.set_pixel(i,j, max(r,new_r),max(g,new_g),max(b,new_b));
	}
	}
}

void colour_in_particle(dla_particle_t& p, unsigned char r, unsigned char g, unsigned char b)
{
	// make the particle at p have colour (r,g,b)
	const int x = p.x;
	const int y = p.y;
	int i,j;
	for (i=x-particle_radius; i<=x+particle_radius; ++i) {
	for (j=y-particle_radius; j<=y+particle_radius; ++j) {
		if (!is_in_image(i,j)) continue;
		double dist = sqrt((double)((x-i)*(x-i)+(y-j)*(y-j)));
		if (dist>particle_radius) continue;

		image.set_pixel(i,j, r,g,b);
	}
	}
}

void save_image_tiled (bitmap_image& full_image, string image_name)
{
    // save image as individual tiles
    int i;
    int j;
    for (i=0; i<tiles_sqrt; ++i) {
    for (j=0; j<tiles_sqrt; ++j) {
        tile.set_all_channels(0,0,0);

        const int64_t x_start = i*tile_x_res;
        const int64_t y_start = j*tile_y_res;
        int64_t x,y;
        for (x=0; x<tile_x_res; ++x) {
        for (y=0; y<tile_y_res; ++y) {
            unsigned char r,g,b;
            full_image.get_pixel(x_start+x,y_start+y, r,g,b);
            tile.set_pixel(x,y, r,g,b);
        }
        }
        tile.save_image(image_name + "_" + to_string(i) + "_" + to_string(j) + ".bmp");
    }
    }
}

pair<double,double> cart_to_polar(int64_t x_int, int64_t y_int)
{
    double x,y,r, theta;
    x = x_int-o_x_dbl;
    y = y_int-o_y_dbl;
    r = sqrt(x*x+y*y);
    theta = atan(y/x);
    pair<double,double> ret {r,theta};
    return ret;
}

int run_model(const int64_t model_iteration)
{
    if (sizeof(void*) != 8)
        {
            cout << "A 64-bit compile/execution is required! Exiting..." << endl;
            return 1;
        }
    assert(particle_radius<med_jump_radius);
    assert(med_jump_radius<long_jump_radius);

	cout << setprecision(2) << fixed;

	list<dla_particle_t> dla_particles;

	// book-keeping:
	// we use the red channel to mark which pixels dla has covered
	// we use the green/blue channel for bookkeeping, see add_particle_to_dla for details
	// the 'true' colours of the image are set at the very end
	// the particles vector stores the order in which the particles landed

	// colour initial disc
	cout << "setting blank image & initial particle... ";
	image.set_all_channels(0,0,0);
	add_particle_to_dla(o_x, o_y);
	dla_particle_t dla_p {(int)o_x,(int)o_y};
	dla_particles.push_back(dla_p);

    dla_radius = particle_radius;
    world_radius_sqr = pow(dla_radius*padding_factor, 2);

	cout << "done" << endl << "particle radius: " << particle_radius << ", origin: " << o_x << "," << o_y << endl;
	if (write_preview_images)
	{
        save_image_tiled(image, "dla_preview");
        cout << "preview iterations: " << preview_iter << endl;
	}

	// file to dump particle positions
	ofstream f;
	f.open("particle_locs_" + to_string((int)particle_radius) + "_" + to_string(model_iteration) + ".csv");
	f << "particle radius is " << particle_radius << ". writing coordinates (x y r theta) relative to origin" << endl;
	f << "n,x,y,r,theta,r^(7/4)/n\n";

	// prepare for main loop
	int64_t iter = 0; // BMs
	int64_t iter_collate = 0;
	int64_t success = 0; // BMs that hit the dla
	int64_t success_collate = 0;
	int64_t collisions = 0; // BMs that hit the dla in a position that was already hit (due to async)
	int64_t collisions_collate = 0;
	int64_t particles = 0; // particles accepted into the dla
	int64_t particles_collate = 0;

	int64_t next_preview_iter = 0;
	int64_t next_status_report = 0;

	const int64_t datapoints_required = max_worker_threads; // do the first 1000 particles single threaded
	const vector<future<list<free_particle_t>>>::size_type n_futures = datapoints_required;
	vector<future<list<free_particle_t>>> futures;
	futures.resize(n_futures); // one receiver slot per thread

	list<free_particle_t> datapoints; // datapoints received from threads are collated to here; a single datapoint is a list of final points of BMs from that thread

	// main loop
	// perform dla until a particle attaches itself in a location that is off-screen
	// in fact; fix the dla and perform some BMs until they hit
	// collate together the hitting locs and add to the dla in batch; repeat
	// note: image is async read safe but not async write safe
	while (true)
	{
		iter_collate = 0;
		success_collate = 0;
		collisions_collate = 0;
		particles_collate = 0;
		world_radius_sqr = pow(dla_radius*padding_factor+1, 2); // +1 for case of particle_radius==0
		BMs_per_thread_per_collate = max((int64_t)1, min((int64_t)max_BMs_per_thread_per_collate, (int64_t)(PI*dla_radius/(300*(1+2.1*particle_radius)))));
        // magic numbers, increase denominator to lower P[collision|success]
        // low BMs_per_thread_per_collate is better accuracy but then the threading will be less efficient, since proportionally more time is spent collating

		int datapoints_requested = 0; // how many we spawned threads for, so far
		int datapoints_recieved = 0; // how many we got back, so far

		// initialize worker threads with a future for each
		for (auto& fu : futures)
		{
			fu = async(launch::async, thread_manager, dla_radius, world_radius_sqr);
			++datapoints_requested;
		}

		// monitor worker threads
		vector<future<list<free_particle_t>>>::iterator fu_it = futures.begin();
		while (datapoints_recieved < datapoints_required)
		{
			future<list<free_particle_t>>& fu = *fu_it;

			if (fu.valid() && fu.wait_for(chrono::milliseconds(2))==future_status::ready)
			{
				// save completed datapoint
				list<free_particle_t> thread_datapoints = fu.get();
				datapoints.splice(datapoints.end(), thread_datapoints);
				++datapoints_recieved;
			}

			++fu_it;
			if (fu_it==futures.end()) fu_it = futures.begin();
		}

		// attempt to add particles from this batch of data into the dla
		// add particles one by one, check for collisions as we go
		// (all particles in this batch were wanting to attach to the same shape dla region)
		bool finished = false;
		for (free_particle_t p : datapoints) {
			++iter;
			++iter_collate;
			if (!p.hit_dla) continue;

			// snap to grid
			int64_t px = p.px();
			int64_t py = p.py();

			// stop if outside dla region
			if (!is_inside_dla_region(px, py)){
				finished = true;
				cout << endl;
				cout.flush();
				cout << "found particle outside max allowed radius at (" << px << "," << py << ")" << endl;
				if (is_in_image(px,py)) image.set_pixel(px, py, 255,0,0);
				if (write_preview_images) save_image_tiled(image, "dla_preview");
				break;
			}

			// attempt to add to dla
			unsigned char r,g,b;
			image.get_pixel(px,py, r,g,b);
			if (g!=127) {
				// don't add; adding this particle would overlap with an already added particle
				++collisions;
				++collisions_collate;
			}
			else {
                // add
				add_particle_to_dla(px, py);
				dla_particle_t dla_p {(int)px,(int)py};
				dla_particles.push_back(dla_p);
				dla_radius = max(dla_radius, p.radial_displacement()+particle_radius);

				pair<double,double> polar = cart_to_polar(px,py);
				double s = pow(polar.first, 1.75) / (double)particles;
				f << particles << "," << px-o_x << "," << py-o_y << "," << polar.first << "," << polar.second << "," << s << endl;
				++particles;
				++particles_collate;
			}

			++success;
			++success_collate;
		}
		datapoints.clear();

		// status report
		if (iter >= next_status_report) {
			cout << "\r" << success << "/" << iter << " particles (radius " << (int)dla_radius
				<< ", P[reject]=" << (double)(iter_collate-success_collate)/(double)iter_collate << "=(" << iter_collate-success_collate << "/" << iter_collate << ")"
				<< ", P[collision|success]=" << (double)collisions_collate/(double)success_collate << "=(" << collisions_collate << "/" << success_collate << ")"
				<< ", BMs_per_collate=" << max_worker_threads*BMs_per_thread_per_collate << ")    ";
			next_status_report = iter + min_BMs_per_status_report;
		}

		if (write_preview_images && (iter>=next_preview_iter || finished)) {
			next_preview_iter = iter + preview_iter;
			save_image_tiled(image, "dla_preview");
		}

		if (iter_collate>100 && success_collate==0) {
			cout << endl;
			cout.flush();
			cout << "Error: could not hit dla" << endl;
			cout << "dla_radius: " << dla_radius << ", world_radius_sqr: " << world_radius_sqr << endl;
			cout << "BMs_per_thread_per_collate: " << BMs_per_thread_per_collate << endl;
			cout << "datapoints_requested: " << datapoints_requested << ", datapoints_recieved: " << datapoints_recieved << endl;
			image.save_image("dla_preview.bmp");
			break;
		}

		if (finished) break;
	}
    // end of main loop
    f.close();

	if (write_image_files)
    {

        // add colours (to new image)
        // red -> yellow -> green -> cyan -> blue -> purple
        // each '->' is 255 steps
        int64_t particles_per_colour = dla_particles.size()/(5*255);
        cout << particles_per_colour << " particles per colour" << endl;
        if (particles_per_colour==0) {
            cout << "not enough particles to shade (" << 5*255 << " needed, have " << dla_particles.size() << ")" << endl;
            save_image_tiled(image, "dla");
            return 0;
        }

        unsigned char r = 255;
        unsigned char g = 0;
        unsigned char b = 0;
        list<dla_particle_t>::size_type i = 0;
        list<dla_particle_t>::iterator i_iter = dla_particles.begin();

        // red -> yellow
        cout << "writing colours...";
        image.set_all_channels(0,0,0); // background
        while (true) {
            dla_particle_t& p = *i_iter;
            colour_in_particle(p, r,g,b);
            ++i;
            ++i_iter;
            if (i%particles_per_colour==0) {
                if (g==255) break;
                ++g;
            }
        }
        // yellow -> green
        while (true) {
            dla_particle_t& p = *i_iter;
            colour_in_particle(p, r,g,b);
            ++i;
            ++i_iter;
            if (i%particles_per_colour==0) {
                if (r==0) break;
                --r;
            }
        }
        // green -> cyan
        while (true) {
            dla_particle_t& p = *i_iter;
            colour_in_particle(p, r,g,b);
            ++i;
            ++i_iter;
            if (i%particles_per_colour==0) {
                if (b==255) break;
                ++b;
            }
        }
        // cyan -> blue
        while (true) {
            dla_particle_t& p = *i_iter;
            colour_in_particle(p, r,g,b);
            ++i;
            ++i_iter;
            if (i%particles_per_colour==0) {
                if (g==0) break;
                --g;
            }
        }
        // blue -> purple
        while (true) {
            if (i_iter==dla_particles.end()) break;
            dla_particle_t& p = *i_iter;
            colour_in_particle(p, r,g,b);
            ++i;
            ++i_iter;
            if (i%particles_per_colour==0) {
                if (r==255) break;
                ++r;
            }
        }

        cout << "done." << endl;

        cout << "writing image files...";
        save_image_tiled(image, "dla");
        cout << "done." << endl;
    }
    else
    {
        cout << "not writing image files to disk" << endl;
    }

    return 0;
}

int main ()
{
    for(int64_t i=0; i<runs; ++i)
    {
        cout << "run " << i+1 << "/" << runs << endl;
        int ret = run_model(i);
        if (ret==1) break;
        cout << endl;
    }
}
