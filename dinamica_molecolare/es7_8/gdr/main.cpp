#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
// #define _USE_MATH_DEFINES
#define LOG(x) cout<<x<<endl;
using namespace std;

class vec {
public:
	double x; double y; double z;

	double mod() {
		return sqrt(x * x + y * y + z * z);
	}
	void uguale(double val) {
		x = val; y = val; z = val;
	}
	void piu(double val) {
		x += val; y += val; z += val;
	}
	void per(double val) {
		x *= val; y *= val; z *= val;
	}
	void stampa() {
		cout << "x: " << x << "y: " << y << "z: " << z << "\n";
	}
};



int main() {
	srand(2);
	const int N = 1e5;
	const int N_bins = 10;

	const double rho = 0.01;
	double L = cbrt(N / rho);
	LOG(L)

	vec r[N], rp;
	for (int i = 0; i < N; ++i) {
		r[i].x = L * rand() / ((double)RAND_MAX + 1.0);
		r[i].y = L * rand() / ((double)RAND_MAX + 1.0);
		r[i].z = L * rand() / ((double)RAND_MAX + 1.0);

		// cout << r[i].x << "\t" << r[i].y << "\t" << r[i].z << "\t" << r[i].mod() << "\n";
	}

	double de_r = 0.5 * L / N_bins;
	double rp_mod, r_k, freq[N_bins], g;

	for (int k = 0; k < N_bins; ++k) {
		freq[k]++;
	}

	for (int p = 0; p < N; ++p) {//particella centrale
		for (int i = 0; i < N; ++i) {//altre N-1 particelle
			rp.x = r[p].x - r[i].x;
			rp.y = r[p].y - r[i].y;
			rp.z = r[p].z - r[i].z;

			rp.x -= L * rint(rp.x / L);
			rp.y -= L * rint(rp.y / L);
			rp.z -= L * rint(rp.z / L);

			rp_mod = rp.mod();
			for (int k = 0; k < N_bins; ++k) {
				r_k = (2 * k + 1) * de_r / 2;
				if (rp_mod >= r_k - de_r / 2 && rp_mod < r_k + de_r / 2) {
					freq[k]++;
					break;
				}
			}
		}
	}

	ofstream gdr; gdr.open("dati.dat");
	for (int k = 0; k < N_bins; ++k) {
		r_k = (2 * k + 1) * de_r / 2;
		freq[k] /= N;
		g = freq[k] / ((pow(r_k + de_r / 2, 3) - pow(r_k - de_r / 2, 3)) * M_PI * 4.0 / 3.0);
		g /= rho;
		gdr << r_k << "\t" << freq[k] << "\t" << g << "\n";
	}
	gdr.close();

	return 0;
}