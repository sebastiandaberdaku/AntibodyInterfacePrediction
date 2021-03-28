/*
 * voxel.cpp
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

/**
 * Implementation of the voxel object.
 */
#include "voxel.h"

/**
 * Default constructor
 */
voxel::voxel() :
		ix(0), iy(0), iz(0) { }
/**
 * Constructor
 */
voxel::voxel(uint16_t i, uint16_t j, uint16_t k) :
		ix(i), iy(j), iz(k) { }
/**
 * Copy constructor.
 */
voxel::voxel(voxel const & vox) :
		ix(vox.ix), iy(vox.iy), iz(vox.iz) { }
/**
 * Copy assignment operator
 */
voxel & voxel::operator=(voxel const & vox) {
	if (this != &vox) {
		this->ix = vox.ix;
		this->iy = vox.iy;
		this->iz = vox.iz;
	}
	return *this;
}
