/*
 * point3D.cpp
 *
 *  Created on: 01/dic/2014
 *      Author: sebastian
 */

#include "point3D.h"
#include "../Atom/atom.h"


point3D::point3D(atom const & atm) : x(atm.x), y(atm.y), z(atm.z) { }
/**
 * Calculates the square of the distance between this point3D and the center of
 * the argument.
 */
float point3D::s_distance(atom const & atm) const {
	float dx = this->x - atm.x;
	float dy = this->y - atm.y;
	float dz = this->z - atm.z;
	return (dx * dx + dy * dy + dz * dz);
};

