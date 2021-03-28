#include "GeometricMoments.h"
#include "Binomial.h"

#ifdef STATIC_GM
int GeometricMoments::maxOrder = -1;

int GeometricMoments::length = 0;
int GeometricMoments::width = 0;
int GeometricMoments::height = 0;

double GeometricMoments::dv = 0;

point3D GeometricMoments::COG = point3D(0, 0, 0);

array1D GeometricMoments::x_samples;
array1D GeometricMoments::y_samples;
array1D GeometricMoments::z_samples;

array2D GeometricMoments::Ir;
array2D GeometricMoments::Is;
array2D GeometricMoments::It;
#endif

GeometricMoments::GeometricMoments() : voxels(NULL) {}

GeometricMoments::GeometricMoments(int maxOrder) : voxels(NULL) {
	moments.resize(maxOrder + 1);
	for (int r = 0; r <= maxOrder; ++r) {
		moments[r].resize(maxOrder - r + 1);
		for (int s = 0; s <= maxOrder - r; ++s) {
			moments[r][s].resize(maxOrder - r - s + 1);
			for (int t = 0; t <= maxOrder - r - s; ++t) {
				moments[r][s][t] = 0;
			}
		}
	}
}

GeometricMoments::GeometricMoments(array3D const & voxels, double scale, point3D const & COG, int maxOrder) {
#ifdef STATIC_GM
#pragma omp critical
	{
		if (this->maxOrder != maxOrder) {
			this->maxOrder = maxOrder;

			dv = scale;

			length = voxels.size();
			width = voxels[0].size();
			height = voxels[0][0].size();

			this->COG = COG;

			x_samples.clear();
			y_samples.clear();
			z_samples.clear();
			computeSamples();

			Ir.clear();
			Is.clear();
			It.clear();
			computeI();
		}
	}
#else
	this->maxOrder = maxOrder;

	dv = scale;

	length = voxels.size();
	width = voxels[0].size();
	height = voxels[0][0].size();

	this->COG = COG;
	computeSamples();
	computeI();
#endif
	this->voxels = &voxels;
	compute();
}

/**
 * Scales and translates the x, y and z values to the
 * new coordinate system centered at the COG of the
 * input 3D function.
 * The values are mapped inside the unit ball.
 */
void GeometricMoments::computeSamples() {
	x_samples.resize(length);
	y_samples.resize(width);
	z_samples.resize(height);

	for (int i = 0; i < length; ++i)
		x_samples[i] = (i - COG.x) * dv;
	for (int j = 0; j < width; ++j)
		y_samples[j] = (j - COG.y) * dv;
	for (int k = 0; k < height; ++k)
		z_samples[k] = (k - COG.z) * dv;
}
/**
 * This method computes the Ir, Is and It data structures.
 */
void GeometricMoments::computeI() {
	Ir.resize(maxOrder + 1);
	Is.resize(maxOrder + 1);
	It.resize(maxOrder + 1);
	for (int ii = 0; ii <= maxOrder; ++ii) {
		Ir[ii].resize(length);
		for (int i = 0; i < length; ++i) {
			Ir[ii][i] = (pow(x_samples[i] + (dv / 2), (ii + 1))
					- pow(x_samples[i] - (dv / 2), (ii + 1))) / (ii + 1);
		}
		Is[ii].resize(width);
		for (int j = 0; j < width; ++j) {
			Is[ii][j] = (pow(y_samples[j] + (dv / 2), (ii + 1))
					- pow(y_samples[j] - (dv / 2), (ii + 1))) / (ii + 1);
		}
		It[ii].resize(height);
		for (int k = 0; k < height; ++k) {
			It[ii][k] = (pow(z_samples[k] + (dv / 2), (ii + 1))
					- pow(z_samples[k] - (dv / 2), (ii + 1))) / (ii + 1);
		}
	}
}
void GeometricMoments::computeI_alt() {
	Ir.resize(maxOrder + 1);
	Is.resize(maxOrder + 1);
	It.resize(maxOrder + 1);

	array2D A;
	array2D B;
	array2D C;

	A.resize(maxOrder + 1);
	B.resize(maxOrder + 1);
	C.resize(maxOrder + 1);

	A[0] = array1D(length, 1);
	B[0] = array1D(width, 1);
	C[0] = array1D(height, 1);

	Ir[0] = array1D(length, dv);
	Is[0] = array1D(width, dv);
	It[0] = array1D(height, dv);

	for (int ii = 1; ii <= maxOrder; ++ii) {
		A[ii].resize(length);
		Ir[ii].resize(length);
		for (int i = 0; i < length; ++i) {
			A[ii][i] = pow(x_samples[i] + (dv / 2), ii) + (x_samples[i] - (dv / 2)) * A[ii - 1][i];
			Ir[ii][i] = (dv / (ii + 1)) * A[ii][i];
		}
		B[ii].resize(width);
		Is[ii].resize(width);
		for (int j = 0; j < width; ++j) {
			B[ii][j] = pow(y_samples[j] + (dv / 2), ii) + (y_samples[j] - (dv / 2)) * B[ii - 1][j];
			Is[ii][j] = (dv / (ii + 1)) * B[ii][j];
		}
		C[ii].resize(height);
		It[ii].resize(height);
		for (int k = 0; k < height; ++k) {
			C[ii][k] = pow(z_samples[k] + (dv / 2), ii) + (z_samples[k] - (dv / 2)) * C[ii - 1][k];
			It[ii][k] = (dv / (ii + 1)) * C[ii][k];
		}
	}
}

void GeometricMoments::computeI_alt2() {
	Ir.resize(maxOrder + 1);
	Is.resize(maxOrder + 1);
	It.resize(maxOrder + 1);

	for (int ii = 0; ii <= maxOrder; ++ii) {
		Ir[ii].resize(length);
		for (int i = 0; i < length; ++i) {
			Ir[ii][i] = 0;
			for (int n = (ii % 2); n <= ii + 1; n += 2) {
				Ir[ii][i] += binomial(ii + 1, n) * pow(x_samples[i], n) * pow(dv / 2, ii + 1 - n);
			}
			Ir[ii][i] *= 2.0 / (double)(ii + 1);
		}
		Is[ii].resize(width);
		for (int j = 0; j < width; ++j) {
			Is[ii][j] = 0;
			for (int n = (ii % 2); n <= ii + 1; n += 2) {
				Is[ii][j] += binomial(ii + 1, n) * pow(y_samples[j], n) * pow(dv / 2, ii + 1 - n);
			}
			Is[ii][j] *= 2.0 / (double)(ii + 1);
		}
		It[ii].resize(height);
		for (int k = 0; k < height; ++k) {
			It[ii][k] = 0;
			for (int n = (ii % 2); n <= ii + 1; n += 2) {
				It[ii][k] += binomial(ii + 1, n) * pow(z_samples[k], n) * pow(dv / 2, ii + 1 - n);
			}
			It[ii][k] *= 2.0 / (double)(ii + 1);
		}
	}
}
void GeometricMoments::computeI_old() {
	Ir.resize(maxOrder + 1);
	Is.resize(maxOrder + 1);
	It.resize(maxOrder + 1);
	for (int ii = 0; ii <= maxOrder; ++ii) {
		Ir[ii].resize(length);
		for (int i = 0; i < length; ++i) {
			Ir[ii][i] = (pow(x_samples[i] + dv, (ii + 1)) - pow(x_samples[i], (ii + 1))) / (ii + 1);
		}
		Is[ii].resize(width);
		for (int j = 0; j < width; ++j) {
			Is[ii][j] = (pow(y_samples[j] + dv, (ii + 1)) - pow(y_samples[j], (ii + 1))) / (ii + 1);
		}
		It[ii].resize(height);
		for (int k = 0; k < height; ++k) {
			It[ii][k] = (pow(z_samples[k] + dv, (ii + 1)) - pow(z_samples[k], (ii + 1))) / (ii + 1);
		}
	}
}
/**
 * This method computes the exact geometric moments of the scaled and
 * translated input 3D function.
 */
void GeometricMoments::compute() {
	array3D Yisk;
	Yisk.resize(length);
	for (int i = 0; i < length; ++i) {
		Yisk[i].resize(maxOrder + 1);
		for (int s = 0; s <= maxOrder; ++s) {
			Yisk[i][s].resize(height);
			for (int k = 0; k < height; ++k) {
				Yisk[i][s][k] = 0;
				for (int j = 0; j < width; ++j) {
					Yisk[i][s][k] += Is[s][j] * (*voxels)[i][j][k];
				}
			}
		}
	}
	array3D Rrsk;
	Rrsk.resize(maxOrder + 1);
	for (int r = 0; r <= maxOrder; ++r) {
		Rrsk[r].resize(maxOrder - r + 1);
		for (int s = 0; s <= maxOrder - r; ++s) {
			Rrsk[r][s].resize(height);
			for (int k = 0; k < height; ++k) {
				Rrsk[r][s][k] = 0;
				for(int i = 0; i < length; ++i) {
					Rrsk[r][s][k] += Yisk[i][s][k] * Ir[r][i];
				}
			}
		}
	}

	moments.resize(maxOrder + 1);
	for (int r = 0; r <= maxOrder; ++r) {
		moments[r].resize(maxOrder - r + 1);
		for (int s = 0; s <= maxOrder - r; ++s){
			moments[r][s].resize(maxOrder - r - s + 1);
			for(int t = 0; t <= maxOrder - r - s; ++t) {
				moments[r][s][t] = 0;
				for (int k = 0; k < height; ++k) {
					moments[r][s][t] += It[t][k] * Rrsk[r][s][k];
				}
#ifdef PRINT_TEST_G
                cout << "Geometric moment[" << r << ",\t" << s << ",\t" << t << "]\t" << moments[r][s][t] << "\n";
#endif
			}
		}
	}
}
/**
 * Directly computes the scaled and translated geometric moments
 * without using any aiding data structure.
 * Written for debugging purposes.
 */
void GeometricMoments::compute_failsafe() {
	moments.resize(maxOrder + 1);
	for (int r = 0; r <= maxOrder; ++r) {
		moments[r].resize(maxOrder - r + 1);
		for (int s = 0; s <= maxOrder - r; ++s) {
			moments[r][s].resize(maxOrder - r - s + 1);
			for (int t = 0; t <= maxOrder - r - s; ++t) {
				moments[r][s][t] = 0;

				for (int i = 0; i < length; ++i) {
					for (int j = 0; j < width; ++j) {
						for (int k = 0; k < height; ++k) {
							moments[r][s][t] += ((pow((i - COG.x + 0.5) * dv, (r + 1)) - pow((i - COG.x - 0.5) * dv, (r + 1)))
									* (pow((j - COG.y + 0.5) * dv, (s + 1)) - pow((j - COG.y - 0.5) * dv, (s + 1)))
									* (pow((k - COG.z + 0.5) * dv, (t + 1)) - pow((k - COG.z - 0.5) * dv, (t + 1)))
									* (*voxels)[i][j][k])
									/ ((r + 1) * (s + 1) * (t + 1));
						}
					}
				}
#ifdef PRINT_TEST_G
				cout << "Geometric moment[" << r << ",\t" << s << ",\t" << t << "]\t" << moments[r][s][t] << "\n";
#endif
			}
		}
	}
}
void GeometricMoments::compute_failsafe_old() {
	moments.resize(maxOrder + 1);
	for (int r = 0; r <= maxOrder; ++r) {
		moments[r].resize(maxOrder - r + 1);
		for (int s = 0; s <= maxOrder - r; ++s) {
			moments[r][s].resize(maxOrder - r - s + 1);
			for (int t = 0; t <= maxOrder - r - s; ++t) {
				moments[r][s][t] = 0;

				for (int i = 0; i < length; ++i) {
					for (int j = 0; j < width; ++j) {
						for (int k = 0; k < height; ++k) {
							moments[r][s][t] += ((pow((i - COG.x + 1) * dv, (r + 1)) - pow((i - COG.x) * dv, (r + 1)))
									* (pow((j - COG.y + 1) * dv, (s + 1)) - pow((j - COG.y) * dv, (s + 1)))
									* (pow((k - COG.z + 1) * dv, (t + 1)) - pow((k - COG.z) * dv, (t + 1)))
									* (*voxels)[i][j][k])
									/ ((r + 1) * (s + 1) * (t + 1));
						}
					}
				}
#ifdef PRINT_TEST_G
				cout << "Geometric moment[" << r << ",\t" << s << ",\t" << t << "]\t" << moments[r][s][t] << "\n";
#endif
			}
		}
	}

}
