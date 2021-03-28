/*
 * voxel_offset.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_VOXEL_OFFSET_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_VOXEL_OFFSET_H_

#include <inttypes.h>
#include <vector>
using namespace std;

struct voxel;
/**
 * Simple struct defining the offset relative to a voxel.
 */
typedef struct voxel_offset {
	int32_t i, j, k;
	voxel_offset();
	voxel_offset(voxel_offset const & s);
	voxel_offset(voxel const & v);
	voxel_offset(int32_t i, int32_t j, int32_t k);
	voxel_offset& operator=(voxel_offset const & s);
	inline const voxel_offset operator-() const {
		return voxel_offset(-(this->i), -(this->j), -(this->k));
	};
	inline voxel_offset & operator+=(voxel_offset const & rhs) {
		this->i += rhs.i;
		this->j += rhs.j;
		this->k += rhs.k;
		return *this;
	};
	inline const voxel_offset operator+(voxel_offset const & rhs) const {
		return voxel_offset(*this) += rhs;
	};
	inline voxel_offset & operator-=(voxel_offset const & rhs) {
		this->i -= rhs.i;
		this->j -= rhs.j;
		this->k -= rhs.k;
		return *this;
	};
	inline const voxel_offset operator-(voxel_offset const & rhs) const {
		return voxel_offset(*this) -= rhs;
	};


	/**
	 * Coordinate transformation for each of the 26 possible neighbours for a single voxel.
	 * By applying one of the following 26 transformations to a voxel's coordinates we
	 * can obtain the coordinates of it's neighbours.
	 */
	static const vector<voxel_offset> neighbours;
} voxel_offset;

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_VOXEL_OFFSET_H_ */
