/*
 * ZernikeMoments.h
 *
 *  Created on: Jul 8, 2014
 *      Author: sebastian
 */

#ifndef ZERNIKEMOMENTS_H_
#define ZERNIKEMOMENTS_H_

#include "../Geometry/point3D.h"
#include "Binomial.h"
#include "Factorial.h"
#include "GeometricMoments.h"
#include <complex>
#include <tuple>
#include <boost/functional/hash.hpp>
#include <unordered_map>

using namespace std;

typedef complex<double> complexValue;

typedef vector<complexValue> complexArray1D; // 3D array of complex type

typedef vector<vector<complexValue> > complexArray2D; // 3D array of complex type

typedef vector<vector<vector<complexValue> > > complexArray3D; // 3D array of complex type

/**
 * Struct representing the complex coefficient of a moment of order (r, s, t)
 */
typedef struct ComplexCoeff {
	ComplexCoeff(int r, int s, int t, complexValue const & value);
	ComplexCoeff(ComplexCoeff const & cc);
	ComplexCoeff();

	int r, s, t;
	complexValue value;
} ComplexCoeff;

typedef vector<vector<vector<vector<ComplexCoeff> > > > ComplexCoeffArray4D;


/**
 * Class representing the Zernike moments
 */
typedef class ZernikeMoments {
public:
	/**
	 * Constructor of the class.
	 * @param order		maximal order to compute the Zernike moments for
	 * @param gm		pointer to the geometrical moments needed for the calculation.
	 */
	ZernikeMoments(int order, GeometricMoments const * gm);
	/**
	 * Default constructor.
	 */
	ZernikeMoments();
	/**
	 * Returns the requested Zernike moment.
	 * Note that: $\Ohm_{n, l}^{-m}(x) = (-1}^m (\Ohm_{n, l}^m)*$
	 */
	inline complexValue getMoment(int n, int l, int m) const {
		if (m >= 0)
			return zernikeMoments[n][l / 2][m];
		else {
			if (m % 2)
				return -conj(zernikeMoments[n][l / 2][abs(m)]);
			else
				return conj(zernikeMoments[n][l / 2][abs(m)]);
		}
	};

	// ---- debug functions/arguments ----

	void reconstruct(complexArray3D & grid, point3D const & COG, double scale);
	void normalizeGridValues(complexArray3D & grid);
	void checkOrthonormality(int n1, int l1, int m1, int n2, int l2, int m2);

private:
	void computeC();
	void computeQ();
	void computeChi();
	void compute();

	complexArray3D zernikeMoments;	// nomen est omen
	static ComplexCoeffArray4D chi;	// coefficients of the geometric moments
	static array3D q; 				// q coefficients (radial polynomial normalization)
	static array2D c;               // c coefficients (harmonic polynomial normalization)
	static int order;        		// := max{n} according to indexing of Zernike polynomials

	GeometricMoments const * gm;

	// ---- debug functions/arguments ----
	void printGrid(complexArray3D& grid);
	double evalMonomialIntegral(int p, int q, int r, int dim);
} ZernikeMoments;



#endif /* ZERNIKEMOMENTS_H_ */
