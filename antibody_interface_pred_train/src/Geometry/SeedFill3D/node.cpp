/*
 * node.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: sebastian
 */


#include "node.h"

node::node() :
		zl(0), zr(0), y(0), direction_y(1) { }

node::node(node const & n) :
		zl(n.zl), zr(n.zr), y(n.y), direction_y(n.direction_y) { }

node::node(int zl, int zr, int y, int direction_y) :
		zl(zl), zr(zr), y(y), direction_y(direction_y) { }

/**
 * Copy assignment operator
 */
node & node::operator=(node const & n) {
	if (this != &n) { // protect against invalid self-assignment
		this->zl = n.zl;
		this->zr = n.zr;
		this->y = n.y;
		this->direction_y = n.direction_y;
	}
	// by convention, always return *this
	return *this;
}
