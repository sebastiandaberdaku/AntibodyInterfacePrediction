/*
 * rigidTransform.h
 *
 *  Created on: 30/lug/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_SVD3_RIGIDTRANSFORM_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_SVD3_RIGIDTRANSFORM_H_

#include "../point3D.h"
#include "svd3.h"

static inline double determinant(double const A[3][3]) {
	return A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2])
			- A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2])
			+ A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
};

/**
 *	This method finds the optimal Rigid/Euclidean transform in 3D space.
 *	It expects as input three 3D points (X[]) to be transformed on top of
 *	three other points (Y[]). It returns R and t: a rotation matrix and
 *	a translation vector which are the optimal/best rotation and translation
 *	to match the points in X on top of the corresponding points in Y.
 */
static inline void rigid_transform_3D(point3D X[3], point3D Y[3],
		double R[3][3], point3D & t) {

	/* centroid of the first three points */
	point3D cX = (X[0] + X[1] + X[2]) / 3;
	/* centroid of the second three points */
	point3D cY = (Y[0] + Y[1] + Y[2]) / 3;
	/*
	 * translated coordinates
	 * dYT is the transposed of dY
	 */
	double dX[3][3] = {
			{X[0].x - cX.x, X[1].x - cX.x, X[2].x - cX.x},
			{X[0].y - cX.y, X[1].y - cX.y, X[2].y - cX.y},
			{X[0].z - cX.z, X[1].z - cX.z, X[2].z - cX.z}
	};

	double dYT[3][3] = {
			{Y[0].x - cY.x, Y[0].y - cY.y, Y[0].z - cY.z},
			{Y[1].x - cY.x, Y[1].y - cY.y, Y[1].z - cY.z},
			{Y[2].x - cY.x, Y[2].y - cY.y, Y[2].z - cY.z}
	};

	double H[3][3];

	matrix_multiply(dX, dYT, H);

	double U[3][3], S[3][3], V[3][3];

	singular_value_decomposition(H, U, S, V);

	double UT[3][3];

	transpose(U, UT);

	matrix_multiply(V, UT, R);

	/*
	 * Special reflection case.
	 *
	 * There’s a special case when finding the rotation matrix that you
	 * have to take care of. Sometimes the SVD will return a ‘reflection’
	 * matrix, which is numerically correct but is actually nonsense in
	 * real life. This is addressed by checking the determinant of R and
	 * seeing if it’s negative. If it is then the 3rd column of R is
	 * multiplied by -1.
	 */
	/* determinant of matrix R */
	double detR = determinant(R);

    if (detR < 0) {
    	R[0][2] *= -1;
    	R[1][2] *= -1;
    	R[2][2] *= -1;
    }
    /*    t = - R x centroid(X) + centroid(Y) */
    t = cY - cX.transform(R);
};


static inline double root_mean_square_deviation(point3D X1[3], point3D Y[3]) {
	point3D err[3];

	err[0] = X1[0] - Y[0];
	err[1] = X1[1] - Y[1];
	err[2] = X1[2] - Y[2];

	double mean_square_error = 0;

	for (int ii = 0; ii < 3; ++ii) {
		mean_square_error += err[ii].s_norm() / 3;
	}

	return sqrt(mean_square_error);
};

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_SVD3_RIGIDTRANSFORM_H_ */
