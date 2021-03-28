/*
 * SurfacePair.h
 *
 *  Created on: Jan 19, 2017
 *      Author: sebastian
 */

#ifndef SURFACEPAIR_H_
#define SURFACEPAIR_H_


#include "../SurfacePatch/CompactPatchDescriptor.h"


class SurfacePair {
public:
	SurfacePair();
	SurfacePair(CompactPatchDescriptor const & p1, CompactPatchDescriptor const & p2);
	SurfacePair(CompactPatchDescriptor const * p1, CompactPatchDescriptor const * p2);
	SurfacePair(SurfacePair const & sp);
	SurfacePair & operator=(SurfacePair const & sp);

	CompactPatchDescriptor const * p1;
	CompactPatchDescriptor const * p2;
	double distance() const;
	double torsion_angle() const;
	double normal1_angle() const;
	double normal2_angle() const;
	bool compatible(SurfacePair const & sp) const;
	SurfacePair inverse_pair() const;

	inline bool operator==(SurfacePair const & rhs) const {
		return (*this->p1 == *rhs.p1) && (*this->p2 == *rhs.p2);
	};
};



#endif /* SURFACEPAIR_H_ */
