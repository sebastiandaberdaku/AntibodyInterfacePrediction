/*
 * voxel.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_VOXEL_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_VOXEL_H_

#include "voxel_offset.h"
#include <inttypes.h>
#include <ostream>
/**
 * Struct defining a voxel data type.
 * A single voxel is univocally identified by it's coordinates.
 */
#if defined TIGHT_PACKING
#pragma pack(push)  /* push current alignment to stack */
#pragma pack(1)     /* set alignment to 1 byte boundary */
#endif
typedef struct voxel {
	uint16_t ix, iy, iz;
	voxel();
	voxel(uint16_t i, uint16_t j, uint16_t k);
	voxel(voxel const & vox);
	voxel & operator=(voxel const & vox);
	inline uint16_t sDistance_to(voxel const & x) const {
		int32_t dx, dy, dz;
		dx = ix - x.ix;
		dy = iy - x.iy;
		dz = iz - x.iz;
		return (dx*dx + dy*dy + dz*dz);
	}
	friend inline std::ostream & operator<<(std::ostream& os, voxel const & v) {
		os << "[" << v.ix << ", " << v.iy << ", " << v.iz << "]";
		return os;
	};

	inline voxel & operator+=(voxel_offset const & rhs) {
		this->ix += rhs.i;
		this->iy += rhs.j;
		this->iz += rhs.k;
		return *this;
	};
	inline const voxel operator+(voxel_offset const & rhs) const {
		return voxel(*this) += rhs;
	};
	inline voxel & operator-=(voxel_offset const & rhs) {
		this->ix -= rhs.i;
		this->iy -= rhs.j;
		this->iz -= rhs.k;
		return *this;
	};
	inline const voxel operator-(voxel_offset const & rhs) const {
		return voxel(*this) -= rhs;
	};
	inline const voxel_offset operator-(voxel const & rhs) const {
		int32_t i = this->ix - rhs.ix;
		int32_t j = this->iy - rhs.iy;
		int32_t k = this->iz - rhs.iz;
		return voxel_offset(i, j, k);
	};
} voxel;
#if defined TIGHT_PACKING
#pragma pack(pop)   /* restore original alignment from stack */
#endif

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_VOXEL_H_ */
