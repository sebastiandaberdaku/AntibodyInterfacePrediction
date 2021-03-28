/*
 * range.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: sebastian
 */


#include "range.h"


range::range() : zl(0), zr(0), y(0) {}

range::range(range const & r) : zl(r.zl), zr(r.zr), y(r.y) {}

range::range(int l, int r, int y) : zl(l), zr(r), y(y) {}

range & range::operator=(range const & r) {
	if (this != &r) { // protect against invalid self-assignment
		this->zl = r.zl;
		this->zr = r.zr;
		this->y = r.y;
	}
	// by convention, always return *this
	return *this;
}
