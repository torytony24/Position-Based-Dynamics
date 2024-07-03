#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

inline float minmod(float x, float y) {
	if (x*y > 0.f) {
		if (x > 0.f) {
			if (x < y) return x;
			else return y;
		}
		else {
			if (x > y) return x;
			else return y;
		}
	}
	else {
		return 0.f;
	}
}

inline double minmod(double x, double y) {
	if (x*y > 0.0) {
		if (x > 0.0) {
			if (x < y) return x;
			else return y;
		}
		else {
			if (x > y) return x;
			else return y;
		}
	}
	else {
		return 0.0;
	}
}

inline void cubicSpline(float f0, float f1, float d0, float d1, float x, float &f, float &d) {

}

// Normalized coordinate
// 0 <= x <= 1
inline void monotonicCubicSpline(
	double f0, double f1,
	double d0, double d1,
	double x, double h,
	double &f, double &d) {

	double delta, a, b;

	delta = f1 - f0;

	d0 *= h;
	d1 *= h;

	// Make monotonic
	if (delta == 0.0) {
		d0 = 0.0; d1 = 0.0;
	}
	else if (delta > 0.0) {
		if (d0 < 0.0) d0 = 0.0;
		else if (d0 > 3.0*delta)
			d0 = 3.0*delta;
		if (d1 < 0.0)	d1 = 0.0;
		else if (d1 > 3.0*delta)
			d1 = 3.0*delta;
	}
	else { // delta < 0.0
		if (d0 > 0.0) d0 = 0.0;
		else if (d0 < 3.0*delta)
			d0 = 3.0*delta;
		if (d1 > 0.0) d1 = 0.0;
		else if (d1 < 3.0*delta)
			d1 = 3.0*delta;
	}

	a = (d0 + d1) - 2.0*delta;
	b = 3.0*delta - (2.0*d0 + d1);

	f = ((a*x + b)*x + d0)*x + f0;
	d = ((3.0*a*x + 2.0*b)*x + d0) / h;
}

// Normalized coordinate
// 0 <= x <= 1
inline void monotonicCubicSpline(
	float f0, float f1,
	float d0, float d1,
	float x, float h,
	float &f, float &d) {

	float delta, a, b;

	delta = f1 - f0;

	d0 *= h;
	d1 *= h;

	// Make monotonic
	if (delta == 0.f) {
		d0 = 0.f; d1 = 0.f;
	}
	else if (delta > 0.f) {
		if (d0 < 0.f) d0 = 0.f;
		else if (d0 > 3.f*delta)
			d0 = 3.f*delta;
		if (d1 < 0.f)	d1 = 0.f;
		else if (d1 > 3.f*delta)
			d1 = 3.f*delta;
	}
	else { // delta < 0.f
		if (d0 > 0.f) d0 = 0.f;
		else if (d0 < 3.f*delta)
			d0 = 3.f*delta;
		if (d1 > 0.f) d1 = 0.f;
		else if (d1 < 3.f*delta)
			d1 = 3.f*delta;
	}

	a = (d0 + d1) - 2.f*delta;
	b = 3.f*delta - (2.f*d0 + d1);

	f = ((a*x + b)*x + d0)*x + f0;
	d = ((3.f*a*x + 2.f*b)*x + d0) / h;
}

// Normalized coordinate
// 0 <= x <= 1
// 0 <= y <= 1
inline void multiCubic2D(
	double f00, double f10, double f01, double f11,
	double dx00, double dx10, double dx01, double dx11,
	double dy00, double dy10, double dy01, double dy11,
	double x, double y,
	double &f, double &dx, double &dy
) {

	double C31, C13, C30, C21, C12, C03, C20, C11, C02, C10, C01, C00;

	C00 = f00;
	C10 = dx00;
	C01 = dy00;

	C20 = 3.f*(f10 - f00) - dx10 - 2.f*dx00;
	C30 = -2.f*(f10 - f00) + dx10 + dx00;

	C02 = 3.f*(f01 - f00) - dy01 - 2.f*dy00;
	C03 = -2.f*(f01 - f00) + dy01 + dy00;

	C21 = 3.f*f11 - 2.f*dx01 - dx11 - C20 - 3.f*(C03 + C02 + C01 + C00);
	C31 = -2.f*f11 + dx01 + dx11 - C30 + 2.f*(C03 + C02 + C01 + C00);

	C12 = 3.f*f11 - 2.f*dy10 - dy11 - C02 - 3.f*(C30 + C20 + C10 + C00);
	C13 = -2.f*f11 + dy10 + dy11 - C03 + 2.f*(C30 + C20 + C10 + C00);

	C11 = dx01 - C13 - C12 - C10;

	f = C31 * x*x*x*y + C13 * x*y*y*y + C30 * x*x*x + C21 * x*x*y + C12 * x*y*y + C03 * y*y*y + C20 * x*x + C11 * x*y + C02 * y*y + C10 * x + C01 * y + C00;
	dx = 3.f*C31*x*x*y + C13 * y*y*y + 3.f*C30*x*x + 2.f*C21*x*y + C12 * y*y + 2.f*C20*x + C11 * y + C10;
	dy = 3.f*C13*y*y*x + C31 * x*x*x + 3.f*C03*y*y + 2.f*C12*y*x + C21 * x*x + 2.f*C02*y + C11 * x + C01;
}

// Normalized coordinate
// 0 <= x <= 1
// 0 <= y <= 1
inline void multiCubic2D(
	float f00, float f10, float f01, float f11,
	float dx00, float dx10, float dx01, float dx11,
	float dy00, float dy10, float dy01, float dy11,
	float x, float y,
	float &f, float &dx, float &dy
) {

	float C31, C13, C30, C21, C12, C03, C20, C11, C02, C10, C01, C00;

	C00 = f00;
	C10 = dx00;
	C01 = dy00;

	C20 = 3.f*(f10 - f00) - dx10 - 2.f*dx00;
	C30 = -2.f*(f10 - f00) + dx10 + dx00;

	C02 = 3.f*(f01 - f00) - dy01 - 2.f*dy00;
	C03 = -2.f*(f01 - f00) + dy01 + dy00;

	C21 = 3.f*f11 - 2.f*dx01 - dx11 - C20 - 3.f*(C03 + C02 + C01 + C00);
	C31 = -2.f*f11 + dx01 + dx11 - C30 + 2.f*(C03 + C02 + C01 + C00);

	C12 = 3.f*f11 - 2.f*dy10 - dy11 - C02 - 3.f*(C30 + C20 + C10 + C00);
	C13 = -2.f*f11 + dy10 + dy11 - C03 + 2.f*(C30 + C20 + C10 + C00);

	C11 = dx01 - C13 - C12 - C10;

	f = C31 * x*x*x*y + C13 * x*y*y*y + C30 * x*x*x + C21 * x*x*y + C12 * x*y*y + C03 * y*y*y + C20 * x*x + C11 * x*y + C02 * y*y + C10 * x + C01 * y + C00;
	dx = 3.f*C31*x*x*y + C13 * y*y*y + 3.f*C30*x*x + 2.f*C21*x*y + C12 * y*y + 2.f*C20*x + C11 * y + C10;
	dy = 3.f*C13*y*y*x + C31 * x*x*x + 3.f*C03*y*y + 2.f*C12*y*x + C21 * x*x + 2.f*C02*y + C11 * x + C01;
}

// Array index
// 0   -> 1   -> 2   -> 3   -> 4   -> 5   -> 6   -> 7
// 000 -> 001 -> 010 -> 011 -> 100 -> 101 -> 110 -> 111
#define I000 0
#define I001 1
#define I010 2
#define I011 3
#define I100 4
#define I101 5
#define I110 6
#define I111 7

// Normalized coordinate
// 0 <= x <= 1
// 0 <= y <= 1
// 0 <= z <= 1
inline void multiCubic3D(
	double *f,
	double *dx,
	double *dy,
	double *dz,
	double x, double y, double z,
	double &F, double &Fx, double &Fy, double &Fz
) {
	double C000, C100, C010, C001;
	double C110, C011, C101, C200, C020, C002;
	double C111, C210, C201, C120, C021, C102, C012, C300, C030, C003;
	double C310, C301, C130, C031, C103, C013, C211, C121, C112;
	double C311, C131, C113;

	double deltax00 = f[I100] - f[I000];
	double deltax01 = f[I101] - f[I001];
	double deltax10 = f[I110] - f[I010];

	double deltay00 = f[I010] - f[I000];
	double deltay01 = f[I011] - f[I001];
	double deltay10 = f[I110] - f[I100];

	double deltaz00 = f[I001] - f[I000];
	double deltaz01 = f[I011] - f[I010];
	double deltaz10 = f[I101] - f[I100];

	C000 = f[I000];
	C100 = dx[I000];
	C010 = dy[I000];
	C001 = dz[I000];

	C310 = dx[I110] - dx[I100] + dx[I010] - dx[I000] - 2.0*(deltax10 - deltax00);
	C210 = 3.0*(deltax10 - deltax00) - 2.0*(dx[I010] - dx[I000]) - (dx[I110] - dx[I100]);

	C301 = dx[I101] - dx[I100] + dx[I001] - dx[I000] - 2.0*(deltax01 - deltax00);
	C201 = 3.0*(deltax01 - deltax00) - 2.0*(dx[I001] - dx[I000]) - (dx[I101] - dx[I100]);

	C130 = dy[I110] - dy[I010] + dy[I100] - dy[I000] - 2.0*(deltay10 - deltay00);
	C120 = 3.0*(deltay10 - deltay00) - 2.0*(dy[I100] - dy[I000]) - (dy[I110] - dy[I010]);

	C031 = dy[I011] - dy[I010] + dy[I001] - dy[I000] - 2.0*(deltay01 - deltay00);
	C021 = 3.0*(deltay01 - deltay00) - 2.0*(dy[I001] - dy[I000]) - (dy[I011] - dy[I010]);

	C103 = dz[I101] - dz[I001] + dz[I100] - dz[I000] - 2.0*(deltaz10 - deltaz00);
	C102 = 3.0*(deltaz10 - deltaz00) - 2.0*(dz[I100] - dz[I000]) - (dz[I101] - dz[I001]);

	C013 = dz[I011] - dz[I001] + dz[I010] - dz[I000] - 2.0*(deltaz01 - deltaz00);
	C012 = 3.0*(deltaz01 - deltaz00) - 2.0*(dz[I010] - dz[I000]) - (dz[I011] - dz[I001]);

	C300 = dx[I100] + dx[I000] - 2.0*deltax00;
	C200 = 3.0*deltax00 - dx[I100] - 2.0*dx[I000];

	C030 = dy[I010] + dy[I000] - 2.0*deltay00;
	C020 = 3.0*deltay00 - dy[I010] - 2.0*dy[I000];

	C003 = dz[I001] + dz[I000] - 2.0*deltaz00;
	C002 = 3.0*deltaz00 - dz[I001] - 2.0*dz[I000];

	C110 = dx[I010] - C100 - C120 - C130;
	C011 = dy[I001] - C010 - C012 - C013;
	C101 = dz[I100] - C001 - C201 - C301;

	double A = f[I100] + dy[I100] + dz[I100] + C011 + C020 + C002 + C120 + C021 + C102 + C012 + C030 + C003 + C130 + C031 + C103 + C013;

	double x0, x1, y0, y1, z0, z1;
	double f111_A = f[I111] - A;

	x0 = dx[I111] - dx[I110] - dx[I101] + dx[I100];
	x1 = dx[I011] - dx[I010] - dx[I001] + dx[I000];
	C311 = x0 + x1 - 2.0*f111_A;
	C211 = 3.0*f111_A - x0 - 2.0*x1;

	y0 = dy[I111] - dy[I110] - dy[I011] + dy[I010];
	y1 = dy[I101] - dy[I100] - dy[I001] + dy[I000];
	C131 = y0 + y1 - 2.0*f111_A;
	C121 = 3.0*f111_A - y0 - 2.0*y1;

	z0 = dz[I111] - dz[I101] - dz[I011] + dz[I001];
	z1 = dz[I110] - dz[I100] - dz[I010] + dz[I000];
	C113 = z0 + z1 - 2.0*f111_A;
	C112 = 3.0*f111_A - z0 - 2.0*z1;

	C111 = x1 + y1 + z1 - 2.0*(f111_A);

	F = C000 + (C001 + (C002 + C003 * z)*z)*z + (C010 + (C011 + (C012 + C013 * z)*z)*z)*y + (C020 + C021 * z)*y*y + (C030 + C031 * z)*y*y*y
		+ (C100 + (C110 + (C120 + C130 * y)*y)*y + (C101 + (C111 + (C121 + C131 * y)*y)*y)*z + (C102 + C112 * y)*z*z + (C103 + C113 * y)*z*z*z)*x
		+ (C200 + C210 * y + (C201 + C211 * y)*z)*x*x
		+ (C300 + C310 * y + (C301 + C311 * y)*z)*x*x*x;

	Fx = C100 + (C110 + (C120 + C130 * y)*y)*y + (C101 + (C111 + (C121 + C131 * y)*y)*y)*z + (C102 + C112 * y)*z*z + (C103 + C113 * y)*z*z*z
		+ 2.0*(C200 + C210 * y + (C201 + C211 * y)*z)*x
		+ 3.0*(C300 + C310 * y + (C301 + C311 * y)*z)*x*x;

	Fy = C010 + (C011 + (C012 + C013 * z)*z)*z + (C110 + (C111 + (C112 + C113 * z)*z)*z)*x + (C210 + C211 * z)*x*x + (C310 + C311 * z)*x*x*x
		+ 2.0*(C020 + C120 * x + (C021 + C121 * x)*z)*y
		+ 3.0*(C030 + C130 * x + (C031 + C131 * x)*z)*y*y;

	Fz = C001 + (C011 + (C021 + C031 * y)*y)*y + (C101 + (C111 + (C121 + C131 * y)*y)*y)*x + (C201 + C211 * y)*x*x + (C301 + C311 * y)*x*x*x
		+ 2.0*(C002 + C102 * x + (C012 + C112 * x)*y)*z
		+ 3.0*(C003 + C103 * x + (C013 + C113 * x)*y)*z*z;
}

// Normalized coordinate
// 0 <= x <= 1
// 0 <= y <= 1
// 0 <= z <= 1
inline void multiCubic3D(
	float *f,
	float *dx,
	float *dy,
	float *dz,
	float x, float y, float z,
	float &F, float &Fx, float &Fy, float &Fz
) {
	float C000, C100, C010, C001;
	float C110, C011, C101, C200, C020, C002;
	float C111, C210, C201, C120, C021, C102, C012, C300, C030, C003;
	float C310, C301, C130, C031, C103, C013, C211, C121, C112;
	float C311, C131, C113;

	float deltax00 = f[I100] - f[I000];
	float deltax01 = f[I101] - f[I001];
	float deltax10 = f[I110] - f[I010];

	float deltay00 = f[I010] - f[I000];
	float deltay01 = f[I011] - f[I001];
	float deltay10 = f[I110] - f[I100];

	float deltaz00 = f[I001] - f[I000];
	float deltaz01 = f[I011] - f[I010];
	float deltaz10 = f[I101] - f[I100];

	C000 = f[I000];
	C100 = dx[I000];
	C010 = dy[I000];
	C001 = dz[I000];

	C310 = dx[I110] - dx[I100] + dx[I010] - dx[I000] - 2.f*(deltax10 - deltax00);
	C210 = 3.f*(deltax10 - deltax00) - 2.f*(dx[I010] - dx[I000]) - (dx[I110] - dx[I100]);

	C301 = dx[I101] - dx[I100] + dx[I001] - dx[I000] - 2.f*(deltax01 - deltax00);
	C201 = 3.f*(deltax01 - deltax00) - 2.f*(dx[I001] - dx[I000]) - (dx[I101] - dx[I100]);

	C130 = dy[I110] - dy[I010] + dy[I100] - dy[I000] - 2.f*(deltay10 - deltay00);
	C120 = 3.f*(deltay10 - deltay00) - 2.f*(dy[I100] - dy[I000]) - (dy[I110] - dy[I010]);

	C031 = dy[I011] - dy[I010] + dy[I001] - dy[I000] - 2.f*(deltay01 - deltay00);
	C021 = 3.f*(deltay01 - deltay00) - 2.f*(dy[I001] - dy[I000]) - (dy[I011] - dy[I010]);

	C103 = dz[I101] - dz[I001] + dz[I100] - dz[I000] - 2.f*(deltaz10 - deltaz00);
	C102 = 3.f*(deltaz10 - deltaz00) - 2.f*(dz[I100] - dz[I000]) - (dz[I101] - dz[I001]);

	C013 = dz[I011] - dz[I001] + dz[I010] - dz[I000] - 2.f*(deltaz01 - deltaz00);
	C012 = 3.f*(deltaz01 - deltaz00) - 2.f*(dz[I010] - dz[I000]) - (dz[I011] - dz[I001]);

	C300 = dx[I100] + dx[I000] - 2.f*deltax00;
	C200 = 3.f*deltax00 - dx[I100] - 2.f*dx[I000];

	C030 = dy[I010] + dy[I000] - 2.f*deltay00;
	C020 = 3.f*deltay00 - dy[I010] - 2.f*dy[I000];

	C003 = dz[I001] + dz[I000] - 2.f*deltaz00;
	C002 = 3.f*deltaz00 - dz[I001] - 2.f*dz[I000];

	C110 = dx[I010] - C100 - C120 - C130;
	C011 = dy[I001] - C010 - C012 - C013;
	C101 = dz[I100] - C001 - C201 - C301;

	float A = f[I100] + dy[I100] + dz[I100] + C011 + C020 + C002 + C120 + C021 + C102 + C012 + C030 + C003 + C130 + C031 + C103 + C013;

	float x0, x1, y0, y1, z0, z1;
	float f111_A = f[I111] - A;

	x0 = dx[I111] - dx[I110] - dx[I101] + dx[I100];
	x1 = dx[I011] - dx[I010] - dx[I001] + dx[I000];
	C311 = x0 + x1 - 2.f*f111_A;
	C211 = 3.f*f111_A - x0 - 2.f*x1;

	y0 = dy[I111] - dy[I110] - dy[I011] + dy[I010];
	y1 = dy[I101] - dy[I100] - dy[I001] + dy[I000];
	C131 = y0 + y1 - 2.f*f111_A;
	C121 = 3.f*f111_A - y0 - 2.f*y1;

	z0 = dz[I111] - dz[I101] - dz[I011] + dz[I001];
	z1 = dz[I110] - dz[I100] - dz[I010] + dz[I000];
	C113 = z0 + z1 - 2.f*f111_A;
	C112 = 3.f*f111_A - z0 - 2.f*z1;

	C111 = x1 + y1 + z1 - 2.f*(f111_A);

	F = C000 + (C001 + (C002 + C003 * z)*z)*z + (C010 + (C011 + (C012 + C013 * z)*z)*z)*y + (C020 + C021 * z)*y*y + (C030 + C031 * z)*y*y*y
		+ (C100 + (C110 + (C120 + C130 * y)*y)*y + (C101 + (C111 + (C121 + C131 * y)*y)*y)*z + (C102 + C112 * y)*z*z + (C103 + C113 * y)*z*z*z)*x
		+ (C200 + C210 * y + (C201 + C211 * y)*z)*x*x
		+ (C300 + C310 * y + (C301 + C311 * y)*z)*x*x*x;

	Fx = C100 + (C110 + (C120 + C130 * y)*y)*y + (C101 + (C111 + (C121 + C131 * y)*y)*y)*z + (C102 + C112 * y)*z*z + (C103 + C113 * y)*z*z*z
		+ 2.f*(C200 + C210 * y + (C201 + C211 * y)*z)*x
		+ 3.f*(C300 + C310 * y + (C301 + C311 * y)*z)*x*x;

	Fy = C010 + (C011 + (C012 + C013 * z)*z)*z + (C110 + (C111 + (C112 + C113 * z)*z)*z)*x + (C210 + C211 * z)*x*x + (C310 + C311 * z)*x*x*x
		+ 2.f*(C020 + C120 * x + (C021 + C121 * x)*z)*y
		+ 3.f*(C030 + C130 * x + (C031 + C131 * x)*z)*y*y;

	Fz = C001 + (C011 + (C021 + C031 * y)*y)*y + (C101 + (C111 + (C121 + C131 * y)*y)*y)*x + (C201 + C211 * y)*x*x + (C301 + C311 * y)*x*x*x
		+ 2.f*(C002 + C102 * x + (C012 + C112 * x)*y)*z
		+ 3.f*(C003 + C103 * x + (C013 + C113 * x)*y)*z*z;
}

#endif // __INTERPOLATION_H__
