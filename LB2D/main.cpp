#include "LBM_2D.h"
#include <ctime>
using namespace std;

void layout_NyNxd() {
	double *f0 = new double[size_ndir];
	double *f1 = new double[size_ndir - size_scalar];
	double *f2 = new double[size_ndir - size_scalar];
	double *rho = new double[size_scalar];
	double *ux = new double[size_scalar];
	double *uy = new double[size_scalar];

	taylor_green(0, rho, ux, uy);
	init_equilibrium(f0, f1, rho, ux, uy);
	clock_t start = clock();
	for (unsigned int n = 1; n <= NSTEPS; n++) {
		bool save = n % 500 == 0;
		//bool save = false;
		stream_collide(f0, f1, f2, rho, ux, uy, save);
		if (save) {
			double prop[4];
			compute_flow_properties(n, rho, ux, uy, prop);
			cout << "t = " << n << endl;
			cout << "energy = " << prop[0] << endl;
			cout << "L2 error in rho = " << prop[1] << endl;
			cout << "L2 error in ux = " << prop[2] << endl;
			cout << "L2 error in uy = " << prop[3] << endl;
		}
		double *temp = f1;
		f1 = f2;
		f2 = temp;
	}

	double runtime = double(clock() - start) / CLOCKS_PER_SEC;
	size_t nodes_updated = NSTEPS * size_t(NX * NY);
	double speed = nodes_updated / (1e6 * runtime);
	double bytesPerGiB = 1024. * 1024. * 1024.;
	double bandwidth = nodes_updated * 2 * ndir * sizeof(double) / (runtime * bytesPerGiB);
	cout << "speed = " << speed << endl;

	delete[]f0;
	delete[]f1;
	delete[]f2;
	delete[]rho;
	delete[]ux;
	delete[]uy;
}

void layout_NyNxd2() {
	double *f0 = new double[size_ndir];
	double *f1 = new double[size_ndir - size_scalar];
	double *rho = new double[size_scalar];
	double *ux = new double[size_scalar];
	double *uy = new double[size_scalar];

	taylor_green(0, rho, ux, uy);
	init_equilibrium(f0, f1, rho, ux, uy);
	clock_t start = clock();
	for (unsigned int n = 1; n <= NSTEPS; n++) {
		bool save = n % 500 == 0;
		//bool save = false;
		stream_collide(f0, f1, rho, ux, uy, save);
		if (save) {
			double prop[4];
			compute_flow_properties(n, rho, ux, uy, prop);
			cout << "t = " << n << endl;
			cout << "energy = " << prop[0] << endl;
			cout << "L2 error in rho = " << prop[1] << endl;
			cout << "L2 error in ux = " << prop[2] << endl;
			cout << "L2 error in uy = " << prop[3] << endl;
		}
	}

	double runtime = double(clock() - start) / CLOCKS_PER_SEC;
	size_t nodes_updated = NSTEPS * size_t(NX * NY);
	double speed = nodes_updated / (1e6 * runtime);
	double bytesPerGiB = 1024. * 1024. * 1024.;
	double bandwidth = nodes_updated * 2 * ndir * sizeof(double) / (runtime * bytesPerGiB);
	cout << "speed = " << speed << endl;

	delete[]f0;
	delete[]f1;
	delete[]rho;
	delete[]ux;
	delete[]uy;
}

int main() {
	layout_NyNxd();
	layout_NyNxd2();
}

