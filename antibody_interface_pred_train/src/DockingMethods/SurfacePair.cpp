/*
 * SurfacePair.cpp
 *
 *  Created on: Jan 19, 2017
 *      Author: sebastian
 */

#include "SurfacePair.h"


SurfacePair::SurfacePair() : p1(NULL), p2(NULL) { };

SurfacePair::SurfacePair(CompactPatchDescriptor const & p1,
		CompactPatchDescriptor const & p2) : p1(&p1), p2(&p2) { };

SurfacePair::SurfacePair(CompactPatchDescriptor const * p1,
		CompactPatchDescriptor const * p2) : p1(p1), p2(p2) { };

SurfacePair::SurfacePair(SurfacePair const & sp) : p1(sp.p1), p2(sp.p2) { };


SurfacePair & SurfacePair::operator=(SurfacePair const & sp) {
	if (this != &sp) {
		this->p1 = sp.p1;
		this->p2 = sp.p2;
	}
	return *this;
};

double SurfacePair::distance() const{
	return (*p1).center.distance((*p2).center);
};
//
//double SurfacePair::torsion_angle() const{
//	return (*p1).patchNormal.angle((*p2).patchNormal);
//};
//
//double SurfacePair::normal1_angle() const{
//	return (*p1).patchNormal.angle((*p2).center - (*p1).center);
//};
//
//double SurfacePair::normal2_angle() const{
//	return (*p2).patchNormal.angle((*p1).center - (*p2).center);
//};

SurfacePair SurfacePair::inverse_pair() const{
	return SurfacePair(this->p2, this->p1);
};
