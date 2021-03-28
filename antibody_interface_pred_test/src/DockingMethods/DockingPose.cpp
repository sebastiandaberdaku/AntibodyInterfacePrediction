/*
 * DockingPose.cpp
 *
 *  Created on: Feb 25, 2015
 *      Author: sebastian
 */
#include "../Geometry/svd3/rigidTransform.h"
#include "../utils/doubleCompare.h"
#include "DockingPose.h"

/**
 * Default constructor.
 */
DockingPose::DockingPose() :
		translation(0, 0, 0) {
}
/**
 * Constructor.
 * @param rotation		the 3x3 rotation matrix
 * @param translation	the translation vector
 */
DockingPose::DockingPose(double const rotation[3][3],
		point3D const & translation) :
		translation(translation) {
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			this->rotation[i][j] = rotation[i][j];
}
/**
 * Copy constructor.
 * @param pose		the DockingPose to copy
 */
DockingPose::DockingPose(DockingPose const & pose) :
		translation(pose.translation) {
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			this->rotation[i][j] = pose.rotation[i][j];
}

/**
 * Constructs the best docking pose for two given SurfaceTriangles
 * and a permutation of their vertexes (ABC if not specified).
 * See the following code for details.
 *
 * @param t1			first SurfaceTriangle
 * @param t2			second SurfaceTriangle
 * @param alignment		the vertexPermutation for the current match
 */
DockingPose::DockingPose(SurfaceTriangle const & t1, SurfaceTriangle const & t2, vertexPermutation alignment) {
	/** vertices of the first triangle */
	point3D T1[3] = { t1.p1->center, t1.p2->center, t1.p3->center };
	/** vertices of the second triangle */
	point3D T2[3];
	switch (alignment) {
	case ABC:
		T2[0] = t2.p1->center; T2[1] = t2.p2->center; T2[2] = t2.p3->center;
		break;
	case ACB:
		T2[0] = t2.p1->center; T2[1] = t2.p3->center; T2[2] = t2.p2->center;
		break;
	case BAC:
		T2[0] = t2.p2->center; T2[1] = t2.p1->center; T2[2] = t2.p3->center;
		break;
	case BCA:
		T2[0] = t2.p2->center; T2[1] = t2.p3->center; T2[2] = t2.p1->center;
		break;
	case CAB:
		T2[0] = t2.p3->center; T2[1] = t2.p1->center; T2[2] = t2.p2->center;
		break;
	case CBA:
		T2[0] = t2.p3->center; T2[1] = t2.p2->center; T2[2] = t2.p1->center;
		break;
	default:
		throw invalid_argument("DockingPose::DockingPose() - Unknown vertex permutation!");
		break;
	}

	rigid_transform_3D(T2, T1, this->rotation, this->translation);
}

/**
 * Copy assignment operator
 */
DockingPose & DockingPose::operator=(DockingPose const & pose) {
	if (this != &pose) {
		this->translation = pose.translation;
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				this->rotation[i][j] = pose.rotation[i][j];
	}
	return *this;
}
ostream & operator<<(ostream & os, DockingPose const & p) {
	os << "Rotation:\t[" << p.rotation[0][0] << ", " << p.rotation[0][1] << ", "
			<< p.rotation[0][2] << "]\n";
	os << "\t\t[" << p.rotation[1][0] << ", " << p.rotation[1][1] << ", "
			<< p.rotation[1][2] << "]\n";
	os << "\t\t[" << p.rotation[2][0] << ", " << p.rotation[2][1] << ", "
			<< p.rotation[2][2] << "]\n";
	os << "Translation:\t" << p.translation << "\n";
	return os;
}

bool DockingPose::operator==(DockingPose const & rhs) const{
	return (this->translation == rhs.translation
			&& doubleCompare(this->rotation[0][0], rhs.rotation[0][0])
			&& doubleCompare(this->rotation[0][1], rhs.rotation[0][1])
			&& doubleCompare(this->rotation[0][2], rhs.rotation[0][2])
			&& doubleCompare(this->rotation[1][0], rhs.rotation[1][0])
			&& doubleCompare(this->rotation[1][1], rhs.rotation[1][1])
			&& doubleCompare(this->rotation[1][2], rhs.rotation[1][2])
			&& doubleCompare(this->rotation[2][0], rhs.rotation[2][0])
			&& doubleCompare(this->rotation[2][1], rhs.rotation[2][1])
			&& doubleCompare(this->rotation[2][2], rhs.rotation[2][2]));
}

bool DockingPose::operator!=(DockingPose const & rhs) {
	return !((*this) == rhs);
}


