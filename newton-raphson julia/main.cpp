#include <iostream>
#include <complex>
#include <vector>
#include <chrono>
#include <iomanip>
#include "bitmap_image.hpp"

using namespace std;
using namespace chrono;

const int iterations = 500;

const int x_res = 2*10236;
const int y_res = 2*7087;
const double x_radius = 2.0;
const double y_radius = x_radius * (double)y_res / (double)x_res;
const double x_incr = 2.0 * x_radius / (double)x_res;
const double y_incr = 2.0 * y_radius / (double)y_res;

const double eps = 1e-9;
const complex<double> minus_one = -1.0;
const complex<double> I = sqrt(minus_one);
const vector<complex<double>> roots {1.0, 0.5*(-1.0+sqrt(3.0)*I), 0.5*(-1.0-sqrt(3.0)*I)};

int found_root (complex<double>& z)
{
	unsigned int i;
	for (i=0; i<roots.size(); ++i) {
		const complex<double>& root = roots[i];
		if (abs(z-root)<eps) return i;
	}
	return -1;
}

inline void iterate (complex<double>& z)
{
	static int i;
	static complex<double> z2;
	static complex<double> z3;
	for (i=0; i<iterations; ++i) {
		z2 = z*z;
		z3 = z2*z;
		z -= (z3-1.0)/(3.0*z2);
	}
}

int main()
{
	cout << setprecision(2) << fixed;
	bitmap_image image {x_res, y_res};

	double x,y;
	unsigned int i,j;

	x = -x_radius;
	for (i=0; i<x_res; ++i) {
		y = -y_radius;
		for (j=0; j<y_res; ++j) {
			complex<double> z {x,y};
			iterate(z);

			int root = found_root(z);
			unsigned char r,g,b;
			switch (root)
			{
			case -1:
				r=0;g=0;b=0;
				break;
			case 0:
				r=255;g=0;b=0;
				break;
			case 1:
				r=0;g=255;b=0;
				break;
			case 2:
				r=0;g=0;b=255;
				break;
			default:
				r=255;g=255;b=255;
				break;
			}
			image.set_pixel(i,j, r,g,b);

			y += y_incr;
		}
		x += x_incr;
		if (i%(x_res/100)==0) cout << "\rprogress: " << (double)i/(double)x_res;
	}
	image.save_image("julia.bmp");

	return 0;
}
