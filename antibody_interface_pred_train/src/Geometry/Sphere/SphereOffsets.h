/*
 * SphereOffsets.h
 *
 *  Created on: Jun 14, 2015
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_SPHERE_SPHEREOFFSETS_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_SPHERE_SPHEREOFFSETS_H_

#include "../voxel_offset.h"
#include "../voxel.h"
#include "../voxelGrid.h"
#include "DrawSphere.h"
#include <list>

using namespace std;

static inline void get_sphere_offsets(int d_rad, list<voxel_offset> & output_sphere_offsets) {
	if (output_sphere_offsets.size() > 0)
		output_sphere_offsets.clear();

	int dim = 2 * d_rad + 1;
	voxelGrid sphere(dim, dim, dim);
	DrawBall(sphere, d_rad, d_rad, d_rad, d_rad);

	for (int ii = 0; ii < dim; ++ii) {
		for (int jj = 0; jj < dim; ++jj) {
			for (int kk = 0; kk < dim; ++kk) {
				if (sphere.getVoxel(ii, jj, kk)) {
					voxel_offset o(ii - d_rad, jj - d_rad, kk - d_rad);
					output_sphere_offsets.push_back(o);
				}
			} // for kk
		} // for jj
	} // for ii

};



#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_SPHERE_SPHEREOFFSETS_H_ */
