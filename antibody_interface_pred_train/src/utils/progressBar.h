/*
 * progressBar.h
 *
 *  Created on: Sep 8, 2015
 *      Author: sebastian
 */

#ifndef PROGRESSBAR_H_
#define PROGRESSBAR_H_

#include <iostream>

using namespace std;

static inline void progressBar(size_t x, size_t n, int f, int w = 50) {
	if ((x % f == 0) || (x == n)) {

		float ratio = x / (float) n;
		int c = ceil(ratio * w);

		cout << setw(3) << (int) (ratio * 100) << "% [";
		for (int i = 0; i < c; ++i)
			cout << "=";
		for (int i = c; i < w; ++i)
			cout << " ";
		cout << "]";
		if (x != n)
			cout << "\r";
		else
			cout << "\n";
		cout << flush;
	}
}

#endif /* PROGRESSBAR_H_ */
