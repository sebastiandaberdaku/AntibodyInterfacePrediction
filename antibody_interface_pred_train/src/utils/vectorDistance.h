/*
 * vectorDistance.h
 *
 *  Created on: 02/nov/2014
 *      Author: sebastian
 */

#ifndef VECTORDISTANCE_H_
#define VECTORDISTANCE_H_

#include <parallel/algorithm>
#include <assert.h>
#include <math.h>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "doubleCompare.h"
using namespace std;


template<typename T>
static inline T KahanSum(vector<T> X) {
    T sum = 0.0;
    T c = 0.0;                  // A running compensation for lost low-order bits.
    for (auto const & x : X) {
        T y = x - c;     // So far, so good: c is zero.
        T t = sum + y;          // Alas, sum is big, y small, so low-order digits of y are lost.
        c = (t - sum) - y; // (t - sum) recovers the high-order part of y; subtracting y recovers -(low part of y)
        sum = t;           // Algebraically, c should always be zero. Beware overly-aggressive optimizing compilers!
        // Next time around, the lost low part will be added to y in a fresh attempt.
    }
    return sum;
}

template<typename T>
static inline T cosine_similarity(vector<T> const & X, vector<T> const & Y) {
	assert(X.size() == Y.size());

	T dot = 0.0, denom_X = 0.0, denom_Y = 0.0;
	for (size_t i = 0; i < X.size(); ++i) {
		dot += X[i] * Y[i];
		denom_X += X[i] * X[i];
		denom_Y += Y[i] * Y[i];
	}
	if (doubleCompare(denom_X, 0.0) || doubleCompare(denom_Y, 0.0))
		return 0.0;
	return dot / (sqrt(denom_X) * sqrt(denom_Y));
//	size_t n = X.size();
//	vector<T> dot(n), sX(n), sY(n);
//
//	transform(X.begin(), X.end(), Y.begin(), dot.begin(), std::multiplies<T>());
//	transform(X.begin(), sX.end(), [](T const & x) {return pow(x, 2.0);});
//	transform(Y.begin(), sY.end(), [](T const & y) {return pow(y, 2.0);});
//	T num = KahanSum(dot), sdenX = KahanSum(sX), sdenY = KahanSum(sY);
//	return num / (sqrt(sdenX) * sqrt(sdenY));
}

/**
 * This method calculates the Pearson product-moment correlation coefficient
 * (sometimes referred to as the PPMCC or PCC or Pearson's r): a measure of
 * the linear correlation (dependence) between two variables X and Y, giving
 * a value between +1 and −1 inclusive, where 1 is total positive correlation,
 * 0 is no correlation, and −1 is total negative correlation.
 * @param X		the first variable
 * @param Y		the second variable
 * @return		the correlation value
 */
//static inline double correlationPearson(vector<double> const & X, vector<double> const & Y) {
//	size_t n = X.size();
//	assert(n == Y.size());
//
//	double Ex = accumulate(X.begin(), X.end(), (double)0.0) / X.size();
//	double Ey = accumulate(Y.begin(), Y.end(), (double)0.0) / Y.size();
//
//	double Sxx = 0, Sxy = 0, Syy = 0;
//	for (size_t ii = 0; ii < n; ++ii) { //Compute the correlation coefficient.
//		double dx = X[ii] - Ex;
//		double dy = Y[ii] - Ey;
//		Sxx += dx * dx;
//		Syy += dy * dy;
//		Sxy += dx * dy;
//	}
//	/*
//	 * TINY_VALUE could be something like 1e-20 and is used to "compensate"
//	 * perfect correlation case (and avoid special verification).
//	 */
//	double TINY_VALUE = 1e-20;
//	return Sxy / (sqrt(Sxx * Syy) + TINY_VALUE);
//};
static inline double correlationPearson(vector<double> const & X, vector<double> const & Y) {
	int n = X.size();

	double ux = 0, uy = 0, vx= 0, vy = 0, wxy = 0;

	for (int ii = 0; ii < n; ++ii) { //Compute the correlation coefficient.
		ux += X[ii];
		uy += Y[ii];

		vx += X[ii] * X[ii];
		vy += Y[ii] * Y[ii];

		wxy += X[ii] * Y[ii];
	}
	/*
	 * TINY_VALUE could be something like 1e-20 and is used to "compensate"
	 * perfect correlation case (and avoid special verification).
	 */
	double TINY_VALUE = 1e-20;
	return (n * wxy - ux * uy)
			/ (sqrt((n * vx - ux * ux) * (n * vy - uy * uy)) + TINY_VALUE);
};
/**
 * This method calculates the squared Euclidean distance between two vectors in a N-dimensional space.
 * @param X		the first vector
 * @param Y		the second vector
 * @return		the squared Euclidean distance between X and Y
 */
static inline double s_distanceEuclidean(vector<double> const & X, vector<double> const & Y) {
	size_t n = X.size();
	assert(n == Y.size());

//	vector<double> temp(n);
//
//	transform(X.begin(), X.end(), Y.begin(), temp.begin(), [](double const & x, double const & y){return pow(x - y, 2.0);});
//	return KahanSum(temp);

	double s_distance = 0; /* squared distance */
	for (size_t i = 0; i < n; ++i) {
		double dxy = X[i] - Y[i];
		s_distance += dxy * dxy;
	}
	return s_distance;
};
/**
 * This method calculates the Euclidean distance between two vectors in a N-dimensional space.
 * In mathematics, the Euclidean distance or Euclidean metric is the "ordinary" distance between
 * two points that one would measure with a ruler, and is given by the Pythagorean formula.
 * @param X		the first vector
 * @param Y		the second vector
 * @return		the Euclidean distance between X and Y
 */
static inline double distanceEuclidean(vector<double> const & X, vector<double> const & Y) {
	return sqrt(s_distanceEuclidean(X, Y));
};
/**
 * This method calculates the Manhattan distance between two vectors in a N-dimensional space.
 * The Manhattan distance between two points is the sum of the absolute differences of their Cartesian
 * coordinates. The Manhattan distance is also known as rectilinear distance, L1 distance, city block
 * distance, taxicab metric, or Manhattan length. The name allude to the grid layout of most streets
 * on the island of Manhattan, which causes the shortest path a car could take between two intersections
 * in the borough to have length equal to the intersections' distance in taxicab geometry.
 * @param X		the first vector
 * @param Y		the second vector
 * @return		the Manhattan distance between X and Y
 */
static inline double distanceManhattan(vector<double> const & X, vector<double> const & Y) {
	size_t n = X.size();
	assert(n == Y.size());
	double distance = 0;
	for (size_t i = 0; i < n; ++i) {
		double dxy = X[i] - Y[i];
		distance += fabs(dxy);
	}
	return distance;
};
template <typename T>
static inline double mean_value(vector<T> const & v) {
	return accumulate(v.begin(), v.end(), (T)0) / (double) v.size();
};
template <typename T>
static inline double sample_standard_deviation(vector<T> const & v) {
	size_t n = v.size();
	if (n < 2)
		throw invalid_argument("At least two elements are required to calculate the standard deviation!");
	double mean = mean_value(v);
	vector<double> diff(n);
	transform(v.begin(), v.end(), diff.begin(), bind2nd(minus<double>(), mean));
	double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	return sqrt(sq_sum / (n - 1));
};
#endif /* VECTORDISTANCE_H_ */
