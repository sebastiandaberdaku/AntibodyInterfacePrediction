/*
 * voxel_offset.cpp
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

/**
 * Implementation of voxel_offset.
 */
#include "voxel.h"
#include "voxel_offset.h"



/**
 * Coordinate transformation for each of the 26 possible neighbours for a single voxel.
 * By applying one of the following 26 transformations to a voxel's coordinates we
 * can obtain the coordinates of it's neighbours.
 */
const vector<voxel_offset> voxel_offset::neighbours = {
		/** Neighbours sharing one face with the current voxel. */
		voxel_offset( 0,  0,  1), voxel_offset( 0,  0, -1), voxel_offset( 0,  1,  0),
		voxel_offset( 0, -1,  0), voxel_offset( 1,  0,  0), voxel_offset(-1,  0,  0),
		/** Neighbours sharing one side with the current voxel. */
		voxel_offset( 1,  1,  0), voxel_offset( 1, -1,  0), voxel_offset(-1,  1,  0),
		voxel_offset(-1, -1,  0), voxel_offset( 1,  0,  1), voxel_offset( 1,  0, -1),
		voxel_offset(-1,  0,  1), voxel_offset(-1,  0, -1), voxel_offset( 0,  1,  1),
		voxel_offset( 0,  1, -1), voxel_offset( 0, -1,  1), voxel_offset( 0, -1, -1),
		/** Neighbours sharing one vertex with the current voxel. */
		voxel_offset( 1,  1,  1), voxel_offset( 1,  1, -1), voxel_offset( 1, -1,  1),
		voxel_offset(-1,  1,  1), voxel_offset( 1, -1, -1), voxel_offset(-1, -1,  1),
		voxel_offset(-1,  1, -1), voxel_offset(-1, -1, -1),
		/** Other neighbours belonging to the 5x5x5 neighbourhood */
		voxel_offset(-2, -2, -2), voxel_offset(-2, -2, -1), voxel_offset(-2, -2,  0),
		voxel_offset(-2, -2,  1), voxel_offset(-2, -2,  2), voxel_offset(-2, -1, -2),
		voxel_offset(-2, -1, -1), voxel_offset(-2, -1,  0), voxel_offset(-2, -1,  1),
		voxel_offset(-2, -1,  2), voxel_offset(-2,  0, -2),	voxel_offset(-2,  0, -1),
		voxel_offset(-2,  0,  0), voxel_offset(-2,  0,  1),	voxel_offset(-2,  0,  2),
		voxel_offset(-2,  1, -2), voxel_offset(-2,  1, -1),	voxel_offset(-2,  1,  0),
		voxel_offset(-2,  1,  1), voxel_offset(-2,  1,  2),	voxel_offset(-2,  2, -2),
		voxel_offset(-2,  2, -1), voxel_offset(-2,  2,  0),	voxel_offset(-2,  2,  1),
		voxel_offset(-2,  2,  2), voxel_offset(-1, -2, -2),	voxel_offset(-1, -2, -1),
		voxel_offset(-1, -2,  0), voxel_offset(-1, -2,  1),	voxel_offset(-1, -2,  2),
		voxel_offset(-1, -1, -2), voxel_offset(-1, -1,  2),	voxel_offset(-1,  0, -2),
		voxel_offset(-1,  0,  2), voxel_offset(-1,  1, -2),	voxel_offset(-1,  1,  2),
		voxel_offset(-1,  2, -2), voxel_offset(-1,  2, -1),	voxel_offset(-1,  2,  0),
		voxel_offset(-1,  2,  1), voxel_offset(-1,  2,  2),	voxel_offset( 0, -2, -2),
		voxel_offset( 0, -2, -1), voxel_offset( 0, -2,  0),	voxel_offset( 0, -2,  1),
		voxel_offset( 0, -2,  2), voxel_offset( 0, -1, -2),	voxel_offset( 0, -1,  2),
		voxel_offset( 0,  0, -2), voxel_offset( 0,  0,  2),	voxel_offset( 0,  1, -2),
		voxel_offset( 0,  1,  2), voxel_offset( 0,  2, -2),	voxel_offset( 0,  2, -1),
		voxel_offset( 0,  2,  0), voxel_offset( 0,  2,  1),	voxel_offset( 0,  2,  2),
		voxel_offset( 1, -2, -2), voxel_offset( 1, -2, -1),	voxel_offset( 1, -2,  0),
		voxel_offset( 1, -2,  1), voxel_offset( 1, -2,  2),	voxel_offset( 1, -1, -2),
		voxel_offset( 1, -1,  2), voxel_offset( 1,  0, -2),	voxel_offset( 1,  0,  2),
		voxel_offset( 1,  1, -2), voxel_offset( 1,  1,  2),	voxel_offset( 1,  2, -2),
		voxel_offset( 1,  2, -1), voxel_offset( 1,  2,  0),	voxel_offset( 1,  2,  1),
		voxel_offset( 1,  2,  2), voxel_offset( 2, -2, -2),	voxel_offset( 2, -2, -1),
		voxel_offset( 2, -2,  0), voxel_offset( 2, -2,  1),	voxel_offset( 2, -2,  2),
		voxel_offset( 2, -1, -2), voxel_offset( 2, -1, -1),	voxel_offset( 2, -1,  0),
		voxel_offset( 2, -1,  1), voxel_offset( 2, -1,  2),	voxel_offset( 2,  0, -2),
		voxel_offset( 2,  0, -1), voxel_offset( 2,  0,  0),	voxel_offset( 2,  0,  1),
		voxel_offset( 2,  0,  2), voxel_offset( 2,  1, -2),	voxel_offset( 2,  1, -1),
		voxel_offset( 2,  1,  0), voxel_offset( 2,  1,  1),	voxel_offset( 2,  1,  2),
		voxel_offset( 2,  2, -2), voxel_offset( 2,  2, -1),	voxel_offset( 2,  2,  0),
		voxel_offset( 2,  2,  1), voxel_offset( 2,  2,  2)
};

/**
 * Default constructor.
 */
voxel_offset::voxel_offset() :
		i(0), j(0), k(0) { }
/**
 * Copy constructor.
 */
voxel_offset::voxel_offset(voxel_offset const & s) :
		i(s.i), j(s.j), k(s.k) { }
/**
 * Offset from voxel.
 */
voxel_offset::voxel_offset(voxel const & v) :
		i(v.ix), j(v.iy), k(v.iz) { }
/**
 * Constructor of the voxel_offset object.
 */
voxel_offset::voxel_offset(int32_t i, int32_t j, int32_t k) :
		i(i), j(j), k(k) { }
/**
 * Copy assignment operator
 */
voxel_offset & voxel_offset::operator=(voxel_offset const & s) {
	if (this != &s) {
		this->i = s.i;
		this->j = s.j;
		this->k = s.k;
	}
	return *this;
}
