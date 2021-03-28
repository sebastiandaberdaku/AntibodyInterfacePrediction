/*
 * Binomial.h
 *
 *  Created on: Jul 8, 2014
 *      Author: sebastian
 */
/**
 * Utility header for the calculation of the binomial coefficient (n choose k).
 * The methods written in this header are numerically safe, in order to calculate
 * high order binomials without numerical overflow.
 */

#ifndef BINOMIAL_H_
#define BINOMIAL_H_

#include <assert.h>
#include <cmath>

using namespace std;

/** Returns the natural logarithm of the binomial "n choose k" */
static inline double logBinomial(int n, int k) {
	assert(k >= 0 && k <= n);

	if (k == 0 || n == k)
		return log(1.0);
	if (k == 1 || n - k == 1)
		return log((double)n);

	double num = 0;
	double den = 0;
	if (n - k < k) {
		for (int ii = k + 1; ii <= n; ++ii) {
			num += log((double) ii);
		}
		for (int ii = 2; ii <= n - k; ++ii) {
			den += log((double) ii);
		}
	} else {
		for (int ii = n - k + 1; ii <= n; ++ii) {
			num += log((double) ii);
		}
		for (int ii = 2; ii <= k; ++ii) {
			den += log((double) ii);
		}
	}
	return (num - den);
};

/** Returns the binomial "n choose k" */
static inline double binomial(int n, int k) {
	return exp(logBinomial(n, k));
};

#endif /* BINOMIAL_H_ */
