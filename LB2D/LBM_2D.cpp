#include "LBM_2D.h"
void inline swap(double &a, double &b) {
	double tmp = a;
	a = b;
	b = tmp;
}

void inline swap(double *f, size_t i, size_t j) {
	double tmp = f[i];
	f[i] = f[j];
	f[j] = tmp;
}

void taylor_green(unsigned int t, unsigned int x, unsigned int y,
	double * r, double * u, double * v) {
	double kx = 2.0 * PI / NX;
	double ky = 2.0 * PI / NY;
	double td = 1.0 / (nu * (kx * kx + ky * ky));
	double X = x + 0.5;
	double Y = y + 0.5;
	double ux = -u_max * sqrt(ky / kx) * cos(kx * X) * sin(ky * Y) * exp(-1.0 * t / td);
	double uy = u_max * sqrt(kx / ky) * sin(kx * X) * cos(ky * Y) * exp(-1.0 * t / td);
	double P = -0.25 * rho0 * u_max * u_max * (
		(ky / kx) * cos(2.0 * kx * X) + (kx / ky) * cos(2.0 * ky * Y)
		) * exp(-2.0 * t / td);
	double rho = rho0 + 3.0 * P;

	*r = rho;
	*u = ux;
	*v = uy;
}

void taylor_green(unsigned int t, double * r, double * u, double * v) {
	for (unsigned int y = 0; y < NY; y++) {
		for (unsigned int x = 0; x < NX; x++) {
			size_t sidx = scalar_index(x, y);
			taylor_green(t, x, y, &r[sidx], &u[sidx], &v[sidx]);
		}
	}
}

void init_equilibrium(double * f, double * r, double * u, double * v) {
	for (unsigned int y = 0; y < NY; y++) {
		for (unsigned int x = 0; x < NX; x++) {
			size_t sidx = scalar_index(x, y);
			double rho = r[sidx];
			double ux = u[sidx];
			double uy = v[sidx];
			double uu = ux * ux + uy * uy;
			for (unsigned int i = 0; i < ndir; i++) {
				double ci_dot_u = dirx[i] * ux + diry[i] * uy;
				f[field_index(x, y, i)] = wi[i] * rho * (1.0 + 3.0 * ci_dot_u + 4.5 * ci_dot_u * ci_dot_u - 1.5 * uu);
			}
		}
	}
}

void init_equilibrium(double * f0, double * f1, double * r, double * u, double * v) {
	for (unsigned int y = 0; y < NY; y++) {
		for (unsigned int x = 0; x < NX; x++) {
			size_t sidx = scalar_index(x, y);
			double rho = r[sidx];
			double ux = u[sidx];
			double uy = v[sidx];
			double uu = ux * ux + uy * uy;
			{
				double ci_dot_u = dirx[0] * ux + diry[0] * uy;
				f0[scalar_index(x, y)] = wi[0] * rho * (1.0 + 3.0 * ci_dot_u + 4.5 * ci_dot_u * ci_dot_u - 1.5 * uu);
			}
			for (unsigned int i = 1; i < ndir; i++) {
				double ci_dot_u = dirx[i] * ux + diry[i] * uy;
				f1[fieldn_index2(x, y, i)] = wi[i] * rho * (1.0 + 3.0 * ci_dot_u + 4.5 * ci_dot_u * ci_dot_u - 1.5 * uu);
			}
		}
	}
}

void compute_rho_u(double * f, double * r, double * u, double * v) {
	for (unsigned int y = 0; y < NY; y++) {
		for (unsigned int x = 0; x < NX; x++) {
			double rho = 0.;
			double ux = 0.;
			double uy = 0.;
			for (unsigned int i = 0; i < ndir; i++) {
				size_t f_idx = field_index(x, y, i);
				rho += f[f_idx];
				ux += dirx[i] * f[f_idx];
				uy += diry[i] * f[f_idx];
			}
			size_t s_idx = scalar_index(x, y);
			r[s_idx] = rho;
			u[s_idx] = ux / rho;
			v[s_idx] = uy / rho;
		}
	}
}

void stream_collide(double * f0, double * f1, double * f2, double * r, double * u, double * v, bool save) {
	for (unsigned int y = 0; y < NY; y++) {
		for (unsigned int x = 0; x < NX; x++) {
			unsigned int xp1 = (x + 1) % NX;
			unsigned int yp1 = (y + 1) % NY;
			unsigned int xm1 = (NX + x - 1) % NX;
			unsigned int ym1 = (NY + y - 1) % NY;
			// direction numbering scheme
			// 6 2 5
			// 3 0 1
			// 7 4 8
			double ft0 = f0[field0_index(x, y)];
			// load populations from adjacent nodes
			double ft1 = f1[fieldn_index2(xm1, y, 1)];
			double ft2 = f1[fieldn_index2(x, ym1, 2)];
			double ft3 = f1[fieldn_index2(xp1, y, 3)];
			double ft4 = f1[fieldn_index2(x, yp1, 4)];
			double ft5 = f1[fieldn_index2(xm1, ym1, 5)];
			double ft6 = f1[fieldn_index2(xp1, ym1, 6)];
			double ft7 = f1[fieldn_index2(xp1, yp1, 7)];
			double ft8 = f1[fieldn_index2(xm1, yp1, 8)];

			// compute moments
			double rho = ft0 + ft1 + ft2 + ft3 + ft4 + ft5 + ft6 + ft7 + ft8;
			double rhoinv = 1.0 / rho;
			double ux = rhoinv*(ft1 + ft5 + ft8 - (ft3 + ft6 + ft7));
			double uy = rhoinv*(ft2 + ft5 + ft6 - (ft4 + ft7 + ft8));

			// only write to memory when needed
			if (save) {
				r[scalar_index(x, y)] = rho;
				u[scalar_index(x, y)] = ux;
				v[scalar_index(x, y)] = uy;
			}

			// now compute and relax to equilibrium
			// temporary variables
			double tw0r = tauinv*w0*rho; // w[0]*rho/tau
			double twsr = tauinv*ws*rho; // w[1-4]*rho/tau
			double twdr = tauinv*wd*rho; // w[5-8]*rho/tau
			double omusq = 1.0 - 1.5*(ux*ux + uy*uy); // 1-(3/2)u.u
			double tux = 3.0*ux;
			double tuy = 3.0*uy;

			f0[field0_index(x, y)] = omtauinv*ft0 + tw0r*(omusq);
			double cidot3u1 = tux;
			f2[fieldn_index2(x, y, 1)] = omtauinv*ft1 + twsr*(omusq + cidot3u1*(1.0 + 0.5*cidot3u1));
			double cidot3u2 = tuy;
			f2[fieldn_index2(x, y, 2)] = omtauinv*ft2 + twsr*(omusq + cidot3u2*(1.0 + 0.5*cidot3u2));
			double cidot3u3 = -tux;
			f2[fieldn_index2(x, y, 3)] = omtauinv*ft3 + twsr*(omusq + cidot3u3*(1.0 + 0.5*cidot3u3));
			double cidot3u4 = -tuy;
			f2[fieldn_index2(x, y, 4)] = omtauinv*ft4 + twsr*(omusq + cidot3u4*(1.0 + 0.5*cidot3u4));
			double cidot3u5 = tux + tuy;
			f2[fieldn_index2(x, y, 5)] = omtauinv*ft5 + twdr*(omusq + cidot3u5*(1.0 + 0.5*cidot3u5));
			double cidot3u6 = -tux + tuy;
			f2[fieldn_index2(x, y, 6)] = omtauinv*ft6 + twdr*(omusq + cidot3u6*(1.0 + 0.5*cidot3u6));
			double cidot3u7 = -tux - tuy;
			f2[fieldn_index2(x, y, 7)] = omtauinv*ft7 + twdr*(omusq + cidot3u7*(1.0 + 0.5*cidot3u7));
			double cidot3u8 = tux - tuy;
			f2[fieldn_index2(x, y, 8)] = omtauinv*ft8 + twdr*(omusq + cidot3u8*(1.0 + 0.5*cidot3u8));
		}
	}
}

void stream(double * f, unsigned int x, unsigned int y) {
	int idx0 = (ndir - 1) * (NX * y + x);
	unsigned int x1 = (x + 1) % NX;
	unsigned int y1 = (y + 1) % NY;
	int idx1 = (ndir - 1) * (NX * y + x1);
	int idx2 = (ndir - 1) * (NX * y1 + x);
	int idx5 = (ndir - 1) * (NX * y1 + x1);
	swap(f, idx0 + 1 - 1, idx1 + 3 - 1);
	swap(f, idx0 + 2 - 1, idx2 + 4 - 1);
	swap(f, idx0 + 5 - 1, idx5 + 7 - 1);
	swap(f, idx1 + 6 - 1, idx2 + 8 - 1);
}

void collide(double * f0, double * f1, unsigned int x, unsigned int y, double *r, double *u, double *v, bool save) {
	double ft0 = f0[field0_index(x, y)];
	// load populations from adjacent nodes
	double *f = f1 + (ndir - 1) * (NX * y + x) - 1;
	double ft1 = f[3];
	double ft2 = f[4];
	double ft3 = f[1];
	double ft4 = f[2];
	double ft5 = f[7];
	double ft6 = f[8];
	double ft7 = f[5];
	double ft8 = f[6];

	// compute moments
	double rho = ft0 + ft1 + ft2 + ft3 + ft4 + ft5 + ft6 + ft7 + ft8;
	double rhoinv = 1.0 / rho;
	double ux = rhoinv*(ft1 + ft5 + ft8 - (ft3 + ft6 + ft7));
	double uy = rhoinv*(ft2 + ft5 + ft6 - (ft4 + ft7 + ft8));

	// only write to memory when needed
	if (save) {
		r[scalar_index(x, y)] = rho;
		u[scalar_index(x, y)] = ux;
		v[scalar_index(x, y)] = uy;
	}

	// now compute and relax to equilibrium
	// temporary variables
	double tw0r = tauinv*w0*rho;              // w[0]*rho/tau
	double twsr = tauinv*ws*rho;              // w[1-4]*rho/tau
	double twdr = tauinv*wd*rho;              // w[5-8]*rho/tau
	double omusq = 1.0 - 1.5*(ux*ux + uy*uy); // 1-(3/2)u.u
	double tux = 3.0*ux;
	double tuy = 3.0*uy;

	f0[field0_index(x, y)] = omtauinv*ft0 + tw0r*(omusq);
	double cidot3u1 = tux;
	f[1] = omtauinv*ft1 + twsr*(omusq + cidot3u1*(1.0 + 0.5*cidot3u1));
	double cidot3u2 = tuy;
	f[2] = omtauinv*ft2 + twsr*(omusq + cidot3u2*(1.0 + 0.5*cidot3u2));
	double cidot3u3 = -tux;
	f[3] = omtauinv*ft3 + twsr*(omusq + cidot3u3*(1.0 + 0.5*cidot3u3));
	double cidot3u4 = -tuy;
	f[4] = omtauinv*ft4 + twsr*(omusq + cidot3u4*(1.0 + 0.5*cidot3u4));
	double cidot3u5 = tux + tuy;
	f[5] = omtauinv*ft5 + twdr*(omusq + cidot3u5*(1.0 + 0.5*cidot3u5));
	double cidot3u6 = -tux + tuy;
	f[6] = omtauinv*ft6 + twdr*(omusq + cidot3u6*(1.0 + 0.5*cidot3u6));
	double cidot3u7 = -tux - tuy;
	f[7] = omtauinv*ft7 + twdr*(omusq + cidot3u7*(1.0 + 0.5*cidot3u7));
	double cidot3u8 = tux - tuy;
	f[8] = omtauinv*ft8 + twdr*(omusq + cidot3u8*(1.0 + 0.5*cidot3u8));
}

void stream_collide(double * f0, double * f1, double *r, double *u, double *v, bool save) {
	for (unsigned int x = 0; x < NX; x++) {
		stream(f1, x, 0);
		stream(f1, x, 1);
	}
	for (unsigned int y = 2; y < NY; y++) {
		for (unsigned int x = 0; x < NX; x++) {
			stream(f1, x, y);
			collide(f0, f1, x, y - 1, r, u, v, save);
		}
	}
	for (unsigned int x = 0; x < NX; x++) {
		collide(f0, f1, x, NY - 1, r, u, v, save);
		collide(f0, f1, x, 0, r, u, v, save);
	}
}

void compute_flow_properties(unsigned int t, double * r, double * u, double * v, double * prop) {
	// prop must point to space for 4 doubles:
	// 0: energy
	// 1: L2 error in rho
	// 2: L2 error in ux
	// 3: L2 error in uy
	double E = 0.0;
	double sumrhoe2 = 0.0;
	double sumuxe2 = 0.0;
	double sumuye2 = 0.0;
	double sumrhoa2 = 0.0;
	double sumuxa2 = 0.0;
	double sumuya2 = 0.0;
	for (unsigned int y = 0; y < NY; ++y) {
		for (unsigned int x = 0; x < NX; ++x) {
			double rho = r[scalar_index(x, y)];
			double ux = u[scalar_index(x, y)];
			double uy = v[scalar_index(x, y)];
			E += rho*(ux*ux + uy*uy);
			double rhoa, uxa, uya;
			taylor_green(t, x, y, &rhoa, &uxa, &uya);
			sumrhoe2 += (rho - rhoa)*(rho - rhoa);
			sumuxe2 += (ux - uxa)*(ux - uxa);
			sumuye2 += (uy - uya)*(uy - uya);
			sumrhoa2 += (rhoa - rho0)*(rhoa - rho0);
			sumuxa2 += uxa*uxa;
			sumuya2 += uya*uya;
		}
	}
	prop[0] = E;
	prop[1] = sqrt(sumrhoe2 / sumrhoa2);
	prop[2] = sqrt(sumuxe2 / sumuxa2);
	prop[3] = sqrt(sumuye2 / sumuya2);
}
