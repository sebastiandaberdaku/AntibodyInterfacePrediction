/*
 * Factorial.h
 *
 *  Created on: Jul 8, 2014
 *      Author: sebastian
 */

/**
 * Utility header for the calculation of the factorial function. All the methods
 * written in this header are numerically safe, allowing to calculate the
 * factorial of large numbers.
 */

#ifndef FACTORIAL_H_
#define FACTORIAL_H_

#include <assert.h>
#include <cmath>

using namespace std;

/**
 * This method implements the Stirling's formula to approximate
 * the value of the factorial function for large numbers.
 */
static inline double stirlingsFactorial(int n) {
		double c = n * log((double) n) + 0.5 * log(2 * M_PI * n) - n;
		return exp(c);
};
/** Gets the natural logarithm of i*(i+1)*...*(j-1)*j */
static inline double logSequenceProduct(int i, int j) {
	assert(i>=1 && j >= i);

	double prod = 0;
	for (int ii = i; ii <= j; ++ii) {
		prod += log((double)ii);
	}

	return prod;
};
/** Gets the natural logarithm of i! */
static inline double logFactorial(int i) {
	if (i == 0 || i == 1)
		return log(1.0);
	return logSequenceProduct(2, i);
};
/** Gets i*(i+1)*...*(j-1)*j */
static inline double sequenceProduct(int i, int j) {
	return exp(logSequenceProduct(i, j));
};
/** Gets i! */
static inline double factorial(int i) {
	return exp(logFactorial(i));
};

#endif /* FACTORIAL_H_ */
