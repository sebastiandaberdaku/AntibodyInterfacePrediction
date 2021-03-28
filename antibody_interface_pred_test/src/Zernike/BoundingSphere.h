/*
 * BoundingSphere.h
 *
 *  Created on: Jul 11, 2015
 *      Author: sebastian
 */

#ifndef BOUNDINGSPHERE_H_
#define BOUNDINGSPHERE_H_

#include <vector>
#include "../Geometry/point3D.h"
#include "GeometricMoments.h"

/**
 * Returns the radius of the tightest bounding sphere centered in COG which contains
 * the whole input 3D image.
 * @param voxels 	3D grid containing the input image
 * @param COG		center of the bounding sphere
 * @return			bounding shpere radius
 */

inline double buondingSphereRadius(array3D const & voxels, point3D const & COG) {
	int length = voxels.size();
	int width = voxels[0].size();
	int height = voxels[0][0].size();
	double max = 0; /* maximum distance from the COG */
	double s_distance = 0;
	double cx = COG.x, cy = COG.y, cz = COG.z;
	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if (fabs(voxels[i][j][k]) > 0) {
					double dx = fabs(i - cx) + 0.5;
					double dy = fabs(j - cy) + 0.5;
					double dz = fabs(k - cz) + 0.5;
					s_distance = dx * dx + dy * dy + dz * dz;
					if (max < s_distance)
						max = s_distance;
				}
			}
		}
	}
	return sqrt(max);
};

inline point3D calculateCOG(array3D const & voxels) {
	int length = voxels.size();
	int width = voxels[0].size();
	int height = voxels[0][0].size();

	size_t count = 0;
	point3D COG(0, 0, 0);
	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if (fabs(voxels[i][j][k]) > 0) {
					++count;
					COG += point3D(i, j, k);
				}
			}
		}
	}
	COG /= count;

	return COG;
};

#endif /* BOUNDINGSPHERE_H_ */
