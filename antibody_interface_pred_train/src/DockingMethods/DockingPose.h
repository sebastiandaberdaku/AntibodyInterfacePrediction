/*
 * DockingPose.h
 *
 *  Created on: 24/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_DOCKINGPOSE_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_DOCKINGPOSE_H_

#include "SurfaceTriangle.h"
#include <boost/functional/hash.hpp>

using namespace std;

/**
 * Struct defining a docking pose between two molecules. A docking pose
 * consists of a rotation matrix and a translation vector.
 * The rotation matrix and the translation vector define the transformation
 * for the current docking pose.
 */
typedef struct DockingPose {
	double rotation[3][3];
	point3D translation;
	/**
	 * Default constructor.
	 */
	DockingPose();

	/**
	 * Constructor.
	 * @param rotation		the 3x3 rotation matrix
	 * @param translation	the translation vector
	 */
	DockingPose(double const rotation[3][3], point3D const & translation);
	/**
	 * Copy constructor.
	 * @param pose		the DockingPose to copy
	 */
	DockingPose(DockingPose const & pose);

	DockingPose(SurfaceTriangle const & t1, SurfaceTriangle const & t2, vertexPermutation alignment = ABC);
	/**
	 * Copy assignment operator
	 */
	DockingPose & operator=(DockingPose const & pose);
	friend ostream & operator<<(ostream & os, DockingPose const & p);
	bool operator==(DockingPose const & rhs) const;
	inline bool operator!=(DockingPose const & rhs);

} DockingPose;

namespace std {
template<>
struct hash<DockingPose> {
	size_t operator()(DockingPose const & p) const {
		size_t result;
		boost::hash_combine(result, p.translation.x);
		boost::hash_combine(result, p.translation.y);
		boost::hash_combine(result, p.translation.z);
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				boost::hash_combine(result, p.rotation[i][j]);
		return result;
	}
};
}

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_DOCKINGPOSE_H_ */
