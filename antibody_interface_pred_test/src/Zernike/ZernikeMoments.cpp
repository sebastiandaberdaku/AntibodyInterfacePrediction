/*
 * ZernikeMoments.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: sebastian
 */

#include "ZernikeMoments.h"
#include "../utils/numerical_limits.h"
#include <assert.h>
#include <iostream>
#include <omp.h>
array2D ZernikeMoments::c;
array3D ZernikeMoments::q;
ComplexCoeffArray4D ZernikeMoments::chi;
int ZernikeMoments::order = -1;

// ---------- implementation of ComplexCoeff struct -------------
/**
 * Copy constructor
 */
ComplexCoeff::ComplexCoeff(ComplexCoeff const & cc) :
		r(cc.r), s(cc.s), t(cc.t), value(cc.value) { }

/**
 * Constructor with scalar args
 */
ComplexCoeff::ComplexCoeff(int r, int s, int t, complexValue const & value) :
		r(r), s(s), t(t), value(value) { }

/**
 * Default constructor
 */
ComplexCoeff::ComplexCoeff() :
		r(0), s(0), t(0), value((0.0, 0.0)) { }


// ---------- implementation of ZernikeMoments class -------------
ZernikeMoments::ZernikeMoments(int order, GeometricMoments const * gm) : gm(gm){
	/**
	 * Computes all coefficients that are input data independent a single time,
	 * as the data structures are declared as static.
	 */
#pragma omp critical
	{
		if (this->order != order) {
			this->order = order;
			c.clear();
			computeC();
			q.clear();
			computeQ();
			chi.clear();
			computeChi();
		}
	}
	compute();
}

ZernikeMoments::ZernikeMoments () : gm(NULL) { }

/**
 * Computes all the normalizing factors $c_l^m$ for harmonic polynomials e
 */
void ZernikeMoments::computeC() {
	/*
	 * indexing:
	 * 	l goes from 0 to n
	 * 	m goes from -l to l, in fact from 0 to l, since c(l,-m) = c(l,m)
	 */
	c.resize(order + 1);
	for (int l = 0; l <= order; ++l) {
		c[l].resize(l + 1);
		for (int m = 0; m <= l; ++m) {
			double log_num = 0.5 * (log(2.0 * l + 1.0) + logFactorial(l + m) + logFactorial(l - m));
			double log_den = logFactorial(l);
			c[l][m] = exp(log_num - log_den);
//			cout << "c["<<l<<"]["<<m<<"]= " << c[l][m]<<endl;
		}
	}
}

/**
 * Computes all coefficients q for the orthonormalization of radial polynomials
 * in Zernike polynomials. Because n and l determine k univocally, it is easier
 * to save q_klnu in q[n][l][nu].
 */
void ZernikeMoments::computeQ() {
	/*
	 * indexing:
	 *	 n  goes from 0 to order
	 *	 l  goes from 0 to n, so that (n - l) is even
	 *	 nu goes from 0 to (n - l) / 2
	 */
	q.resize(order + 1);            // there is order + 1 n's
	for (int n = 0; n <= order; ++n) {
		q[n].resize(n / 2 + 1);      // there is floor(n/2) + 1 l's
		for (int l = n % 2; l <= n; l += 2) {
			int k = (n - l) / 2;
			q[n][l / 2].resize(k + 1);   // there is k+1 nu's
			for (int nu = 0; nu <= k; ++nu) {
				// numerator
				double log_num = logBinomial(2 * k, k) + logBinomial(k, nu)
						+ logBinomial(2 * (k + l + nu) + 1, 2 * k)
						+ 0.5 * log(2.0 * l + 4.0 * k + 3.0);
				// denominator
				double log_den = 2 * k * log(2.0) + logBinomial(k + l + nu, k)
						+ 0.5 * log(3.0);

				q[n][l / 2][nu] = exp(log_num - log_den);
				if ((k + nu) % 2) {
					q[n][l / 2][nu] *= -1;
				}
			}
		}
	}
}

/**
 * Computes the coefficients of geometrical moments in linear combinations
 * yielding the Zernike moments for each applicable [n,l,m] for n<=order.
 * For each such combination the coefficients are stored with according
 * geometrical moment order (see ComplexCoeff).
 */
void ZernikeMoments::computeChi() {
	chi.resize(order + 1);
	for (int n = 0; n <= order; ++n) {
		chi[n].resize(n / 2 + 1); // there is floor(n/2) + 1 l's
		for (int l = n % 2; l <= n; l += 2) {
			chi[n][l / 2].resize(l + 1); // m goes from 0 to l
			for (int m = 0; m <= l; ++m) {
				double w = c[l][m] / pow(2.0, (double)m);
				int k = (n - l) / 2; // 2k = n - l
				for (int nu = 0; nu <= k; ++nu) {
					double wnu = w * q[n][l / 2][nu]; // q_{kl}^nu
					for (int alpha = 0; alpha <= nu; ++alpha) {
						double wnuA = wnu * binomial(nu, alpha);
						for (int beta = 0; beta <= nu - alpha; ++beta) {
							double wnuAB = wnuA * binomial(nu - alpha, beta);
							for (int u = 0; u <= m; ++u) {
								double wnuABu = wnuAB * binomial(m, u);
								for (int mu = 0; mu <= (l - m) / 2; ++mu) {
									double wnuABuMu = wnuABu
											* binomial(l, mu)
											* binomial(l - mu, m + mu)
											/ pow(2.0, 2.0 * mu);
									for (int eta = 0; eta <= mu; ++eta) {
										// the absolute value of the coefficient
										double wnuABuMuEta = wnuABuMu * binomial(mu, eta);
										// the sign
										if ((m - u + mu) % 2)
											wnuABuMuEta *= -1;
										// * i^p
										int rest = u % 4;
										complexValue chi_nlm;
										switch (rest) {
											case 0: chi_nlm = complexValue(wnuABuMuEta, 0.0); break;
											case 1: chi_nlm = complexValue(0.0, wnuABuMuEta); break;
											case 2: chi_nlm = complexValue(-wnuABuMuEta, 0.0); break;
											case 3: chi_nlm = complexValue(0.0, -wnuABuMuEta); break;
										}
										// determination of the order of the according moment
										int r = 2 * (eta + alpha) + u;
										int s = 2 * (mu - eta + beta) + m - u;
										int t = 2 * (nu - alpha - beta - mu) + l - m;
										ComplexCoeff cc(r, s, t, chi_nlm);
										chi[n][l / 2][m].push_back(cc);
//										cout << "chi["<<n<<"]["<<l <<"]["<<m<<"]["<<r<<"]["<<s<<"]["<<t<<"]=" << cc.value <<endl;
									} // eta
								} // mu
							} // u
						} // beta
					} // alpha
				} // nu
			} // m
		} // l
	} // n
}


/**
 * Computes the Zernike moments. This computation is data dependent
 * and has to be performed for each new object and/or transformation.
 */
void ZernikeMoments::compute() {
	// geometric moments have to be computed first
	assert(order > 0);
	assert(gm != NULL);
	// the maximal order of the computed geometric moments must be
	// at least equal to the maximal order of the Zernike moments
	assert(gm->maxOrder >= order);
	/*
	 * indexing:
	 *	 n goes 0..order
	 *	 l goes 0..n so that n-l is even
	 *	 m goes -l..l
	 */
	zernikeMoments.resize(order + 1);
	for (int n = 0; n <= order; ++n) {
		zernikeMoments[n].resize(n / 2 + 1);
		for (int l = n % 2; l <= n; l += 2) {
			zernikeMoments[n][l / 2].resize(l + 1);
			for (int m = 0; m <= l; ++m) {
				// Zernike moment of according indices [nlm]
				zernikeMoments[n][l / 2][m] = complexValue(0.0, 0.0);
				for (auto const & chi_nlm : chi[n][l / 2][m]) {
					zernikeMoments[n][l / 2][m] += conj(chi_nlm.value) * (gm->getMoment(chi_nlm.r, chi_nlm.s, chi_nlm.t));
				}
				zernikeMoments[n][l / 2][m] *= 3.0 / (4.0 * M_PI);
#ifdef PRINT_TEST_Z
                cout << "Zernike moment[" << n << ",\t" << l << ",\t" << m << "]\t" << zernikeMoments[n][l / 2][m] << "\n";
#endif
			}
		}
	}
}

/**
 * This method reconstructs the input 3D function from its complex valued Zernike moments.
 *
 * @param grid		grid that will contain the reconstructed 3D function encoded
 * 					in its complex valued Zernike moments
 * @param COG		center of gravity
 * @param scale		scaling factor to map into the unit ball
 * @param minN		min value for n freq index
 * @param maxN		max value for n freq index
 * @param minL		min value for l freq index
 * @param maxL		max value for l freq index
 */
void ZernikeMoments::reconstruct(complexArray3D& grid, point3D const & COG, double scale) {
	int length = grid.size();
	int width = grid[0].size();
	int height = grid[0][0].size();
	cout << "Reconstructing 3D image inside a "
			<< length << "x" << width << "x" << height << " grid.\n";
	cout << "Scale factor: " << scale << "\n";
	
	double dx = scale, dy = scale, dz = scale;
	
	array1D x_samples(length);
	array1D y_samples(width);
	array1D z_samples(height);

	for (int i = 0; i < length; ++i)
		x_samples[i] = (i - COG.x) * dx;
	for (int j = 0; j < width; ++j)
		y_samples[j] = (j - COG.y) * dy;
	for (int k = 0; k < height; ++k)
		z_samples[k] = (k - COG.z) * dz;

#pragma omp parallel for
	for (int i = 0; i < length; ++i) {
#pragma omp critical
		{
			int id = omp_get_thread_num();
			cout << "Thread["<< id  <<"]: Reconstructing layer: " << i << endl;
		}
        for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
                // the origin is in the middle of the grid, all voxels are
                // projected into the unit ball
			    //current point
				point3D point(x_samples[i], y_samples[j], z_samples[k]);
				if (point.s_norm() > 1.0)
					continue;

				// function value at point
				complexValue fVal(0.0, 0.0);
				for (int n = 0; n <= order; ++n) {
					int maxK = n / 2;
					for (int kk = 0; kk <= maxK; ++kk) {
						for (int nu = 0; nu <= kk; ++nu) {
							int l = n - 2 * kk;
							for (int m = -l; m <= l; ++m) {
								// zernike polynomial evaluated at point
								complexValue zp(0.0, 0.0);

								int nCoeffs = chi[n][l / 2][abs(m)].size();
								for (int ii = 0; ii < nCoeffs; ++ii) {
									ComplexCoeff cc = chi[n][l / 2][abs(m)][ii];
									complexValue cvalue = cc.value;
									double p = cc.r, q = cc.s, r = cc.t;
									// conjugate if m negative
									if (m < 0) {
										cvalue = conj(cvalue);
										// take care of the sign
										if (m % 2)
											cvalue *= -1.0;
									}

									zp += cvalue * pow(point.x, p) * pow(point.y, q) * pow(point.z, r);
								}
								fVal += zp * getMoment(n, l, m);
							}
						}
					}
				}
				grid[i][j][k] = fVal;
			}
		}
	}
	normalizeGridValues(grid);
}

/**
 *
 */
void ZernikeMoments::normalizeGridValues(complexArray3D& grid) {
	int length = grid.size();
	int width = grid[0].size();
	int height = grid[0][0].size();
	cout << "Normalizing 3D complex function image inside a "
			<< length << "x" << width << "x" << height << " grid.\n";
	double max = min_double;
	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if(norm(grid[i][j][k]) > max)
					max = norm(grid[i][j][k]);
			}
		}
	}

	assert(max > 0);
	cout << "\nMaximal value in grid: " << max << "\n";
	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				grid[i][j][k] /= sqrt(max);
			}
		}
	}
}

/**
 * Prints the given compex valued grid for debugging.
 * @param grid		complex valued 3D grid to print
 */
void ZernikeMoments::printGrid(complexArray3D& grid) {
	int xD = grid.size();
	int yD = grid[0].size();
	int zD = grid[0][0].size();

	cout.setf(ios_base::scientific, ios_base::floatfield);

	double max = 0;
	for (int k = 0; k < zD; ++k) {
		for (int j = 0; j < yD; ++j) {
			for (int i = 0; i < xD; ++i) {
				if (fabs(grid[i][j][k].real()) > max) {
					max = fabs(grid[i][j][k].real());
				}
			}
		}
	}

	for (int k = 0; k < zD; ++k) {
		cout << k << ". layer:\n";
		for (int j = 0; j < yD; ++j) {
			for (int i = 0; i < xD; ++i) {
				//cout << grid[i][j][k].real () / max << "\t";
				cout << grid[i][j][k].real() << "\t";
			}
			cout << "\n";
		}
		cout << "\n";
	}
	cout.setf(ios_base::fmtflags(0), ios_base::floatfield);
}
/**
 * Calculates the inner product between Zernike moments
 * of order [n1, l1, m1] and [n2, l2, m2].
 */
void ZernikeMoments::checkOrthonormality(int n1, int l1, int m1,
		int n2, int l2,	int m2) {
	int dim = 64;

	// the total sum of the scalar product
	complexValue sum(0.0, 0.0);

	int nCoeffs1 = (int) chi[n1][l1 / 2][m1].size();
	int nCoeffs2 = (int) chi[n2][l2 / 2][m2].size();

	for (int i = 0; i < nCoeffs1; ++i) {
		ComplexCoeff cc1 = chi[n1][l1 / 2][m1][i];
		for (int j = 0; j < nCoeffs2; ++j) {
			ComplexCoeff cc2 = chi[n2][l2 / 2][m2][j];

			int r = cc1.r + cc2.r;
			int s = cc1.s + cc2.s;
			int t = cc1.t + cc2.t;

			sum += cc1.value * conj(cc2.value) * evalMonomialIntegral(r, s, t, dim);
		}
	}

	cout << "\nInner product of [" << n1 << "," << l1 << "," << m1 << "]";
	cout << " and [" << n2 << "," << l2 << "," << m2 << "]: ";
	cout << sum << "\n\n";
}


/**
 * Evaluates the integral of a monomial x^p*y^q*z^r within the unit ball
 * Attention : a very stupid implementation, thus it's accordingly very slow
 */
double ZernikeMoments::evalMonomialIntegral(int p, int q, int r, int dim) {
	double radius = (dim - 1) / 2.0;
	double scale = pow(1 / radius, 3.0);
	double center = (dim - 1) / 2.0;

	double result = 0.0;
	point3D point;

	for (int x = 0; x < dim; ++x) {
		point.x = (x - center) / radius;
		for (int y = 0; y < dim; ++y) {
			point.y = (y - center) / radius;
			for (int z = 0; z < dim; ++z) {
				point.z = (z - center) / radius;

				if (point.s_norm() > 1) {
					continue;
				}

				result += pow(point.x, (double) p) * pow(point.y, (double) q)
						* pow(point.z, (double) r);
			}
		}
	}

	result *= (3.0 / (4.0 * M_PI)) * scale;
	return result;
}
