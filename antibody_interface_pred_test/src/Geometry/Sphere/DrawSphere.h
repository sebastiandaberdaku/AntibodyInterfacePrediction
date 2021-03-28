/*
 * DrawSphere.h
 *
 *  Created on: 02/feb/2015
 *      Author: sebastian
 */

#ifndef SPHERE_DRAWSPHERE_H_
#define SPHERE_DRAWSPHERE_H_

#include "../voxelGrid.h"

/**
 * Implementation of the Midpoint Circle Algorithm which determines the points needed for
 * drawing a circle on a given plane (parallel to the XY plane). This method gets an initial
 * error value as an additional parameter, in order to compensate for the rounding error of
 * the center coordinates and radius.
 *
 * If no initial error is propagated the sphere would raster to a cylinder near its equator,
 * as the same circle is drawn several times, where in reality the arc should vary slightly
 * even if the x=0 or y=0 points would stay the same.
 *
 * Using the accumulated error from each point as an input to the next circle solves the issue.
 */
inline static void DrawCircleXY(voxelGrid & grid, int x0, int y0, int z, int radius, int error0) {
	int y = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (y >= x) {
		grid.setVoxel(x0 + x, y0 + y, z);
		grid.setVoxel(x0 + x, y0 - y, z);
		grid.setVoxel(x0 - x, y0 + y, z);
		grid.setVoxel(x0 - x, y0 - y, z);
		grid.setVoxel(x0 + y, y0 + x, z);
		grid.setVoxel(x0 + y, y0 - x, z);
		grid.setVoxel(x0 - y, y0 + x, z);
		grid.setVoxel(x0 - y, y0 - x, z);

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--y;
			radiusError += 2 * (x - y) + 1;
		}
	}
};
template<typename T>
inline static void DrawCircleXY(vector<vector<vector<T>>> & grid, T value, int x0, int y0, int z, int radius, int error0) {
	int y = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (y >= x) {
		grid[x0 + x][y0 + y][z] = value;
		grid[x0 + x][y0 - y][z] = value;
		grid[x0 - x][y0 + y][z] = value;
		grid[x0 - x][y0 - y][z] = value;
		grid[x0 + y][y0 + x][z] = value;
		grid[x0 + y][y0 - x][z] = value;
		grid[x0 - y][y0 + x][z] = value;
		grid[x0 - y][y0 - x][z] = value;

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--y;
			radiusError += 2 * (x - y) + 1;
		}
	}
};


inline static void DrawDiskXY(voxelGrid & grid, int x0, int y0, int z, int radius, int error0) {
	int y = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (y >= x) {
		for (int i = y; i >= x; --i) {
			grid.setVoxel(x0 + x, y0 + i, z);
			grid.setVoxel(x0 - x, y0 + i, z);
			grid.setVoxel(x0 - i, y0 + x, z);
			grid.setVoxel(x0 + i, y0 + x, z);

			grid.setVoxel(x0 + x, y0 - i, z);
			grid.setVoxel(x0 - x, y0 - i, z);
			grid.setVoxel(x0 - i, y0 - x, z);
			grid.setVoxel(x0 + i, y0 - x, z);
		}

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--y;
			radiusError += 2 * (x - y) + 1;
		}
	}
};

/**
 * Implementation of the Midpoint Circle Algorithm which determines the points needed for
 * drawing a circle on a given plane (parallel to the YZ plane).
 */
inline static void DrawCircleYZ(voxelGrid & grid, int x, int y0, int z0, int radius, int error0) {
	int z = radius, y = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= y) {
		grid.setVoxel(x, y0 + y, z0 + z);
		grid.setVoxel(x, y0 + y, z0 - z);
		grid.setVoxel(x, y0 - y, z0 + z);
		grid.setVoxel(x, y0 - y, z0 - z);
		grid.setVoxel(x, y0 + z, z0 + y);
		grid.setVoxel(x, y0 + z, z0 - y);
		grid.setVoxel(x, y0 - z, z0 + y);
		grid.setVoxel(x, y0 - z, z0 - y);

		++y;
		if (radiusError < 0) {
			radiusError += 2 * y + 1;
		} else {
			--z;
			radiusError += 2 * (y - z) + 1;
		}
	}
};
template<typename T>
inline static void DrawCircleYZ(vector<vector<vector<T>>> & grid, T value, int x, int y0, int z0, int radius, int error0) {
	int z = radius, y = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= y) {
		grid[x][y0 + y][z0 + z] = value;
		grid[x][y0 + y][z0 - z] = value;
		grid[x][y0 - y][z0 + z] = value;
		grid[x][y0 - y][z0 - z] = value;
		grid[x][y0 + z][z0 + y] = value;
		grid[x][y0 + z][z0 - y] = value;
		grid[x][y0 - z][z0 + y] = value;
		grid[x][y0 - z][z0 - y] = value;

		++y;
		if (radiusError < 0) {
			radiusError += 2 * y + 1;
		} else {
			--z;
			radiusError += 2 * (y - z) + 1;
		}
	}
};
/**
 * Implementation of the Midpoint Circle Algorithm which determines the points needed for
 * drawing a circle on a given plane (parallel to the YZ plane).
 */
inline static void DrawDiskYZ(voxelGrid & grid, int x, int y0, int z0, int radius, int error0) {
	int z = radius, y = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= y) {
		for (int i = z; i >= y; --i) {
			grid.setVoxel(x, y0 + i, z0 + y);
			grid.setVoxel(x, y0 + i, z0 - y);
			grid.setVoxel(x, y0 - i, z0 + y);
			grid.setVoxel(x, y0 - i, z0 - y);

			grid.setVoxel(x, y0 + y, z0 + i);
			grid.setVoxel(x, y0 + y, z0 - i);
			grid.setVoxel(x, y0 - y, z0 + i);
			grid.setVoxel(x, y0 - y, z0 - i);
		}
		++y;
		if (radiusError < 0) {
			radiusError += 2 * y + 1;
		} else {
			--z;
			radiusError += 2 * (y - z) + 1;
		}
	}
};



/**
 * Implementation of the Midpoint Circle Algorithm which determines the points needed for
 * drawing a circle on a given plane (parallel to the XZ plane).
 */
inline static void DrawCircleXZ(voxelGrid & grid, int x0, int y, int z0, int radius, int error0) {
	int z = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= x) {
		grid.setVoxel(x0 + x, y, z0 + z);
		grid.setVoxel(x0 + x, y, z0 - z);
		grid.setVoxel(x0 - x, y, z0 + z);
		grid.setVoxel(x0 - x, y, z0 - z);
		grid.setVoxel(x0 + z, y, z0 + x);
		grid.setVoxel(x0 + z, y, z0 - x);
		grid.setVoxel(x0 - z, y, z0 + x);
		grid.setVoxel(x0 - z, y, z0 - x);

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--z;
			radiusError += 2 * (x - z) + 1;
		}
	}
};
template<typename T>
inline static void DrawCircleXZ(vector<vector<vector<T>>> & grid, T value, int x0, int y, int z0, int radius, int error0) {
	int z = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= x) {
		grid[x0 + x][y][z0 + z] = value;
		grid[x0 + x][y][z0 - z] = value;
		grid[x0 - x][y][z0 + z] = value;
		grid[x0 - x][y][z0 - z] = value;
		grid[x0 + z][y][z0 + x] = value;
		grid[x0 + z][y][z0 - x] = value;
		grid[x0 - z][y][z0 + x] = value;
		grid[x0 - z][y][z0 - x] = value;

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--z;
			radiusError += 2 * (x - z) + 1;
		}
	}
};
inline static void DrawDiskXZ(voxelGrid & grid, int x0, int y, int z0, int radius, int error0) {
	int z = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= x) {
		for (int i = z0 - z; i <= z0 + z; ++i) {
			grid.setVoxel(x0 + x, y, i);
			grid.setVoxel(x0 - x, y, i);
		}
		for (int i = z0 - x; i <= z0 + x; ++i) {
			grid.setVoxel(x0 + z, y, i);
			grid.setVoxel(x0 - z, y, i);
		}

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--z;
			radiusError += 2 * (x - z) + 1;
		}
	}
};
/**
 * Adapted Midpoint Circle Algorithm for drawing spheres in 3D.
 * @param grid		target voxelGrid
 * @param x0		x coordinate of the center of the sphere
 * @param y0		y coordinate of the center of the sphere
 * @param z0		z coordinate of the center of the sphere
 * @param radius	radius of the sphere
 */

inline static void DrawSphere(voxelGrid & grid, int x0, int y0, int z0, int radius) {
	int r = radius, t = 0;
	int radiusError = 1 - r;

	while (r >= t) {
		// pass in base point (x0,y0,z0), this algorithm's t as the current plane,
		// this algorithm's r as the radius, and pass along radius error.
		DrawCircleXY(grid, x0, y0, z0 - t, r, radiusError);
		DrawCircleXY(grid, x0, y0, z0 + t, r, radiusError);
		DrawCircleYZ(grid, x0 - t, y0, z0, r, radiusError);
		DrawCircleYZ(grid, x0 + t, y0, z0, r, radiusError);
		DrawCircleXZ(grid, x0, y0 + t, z0, r, radiusError);
		DrawCircleXZ(grid, x0, y0 - t, z0, r, radiusError);
		++t;
		if (radiusError < 0) {
			radiusError += 2 * t + 1;
		} else {
			r--;
			radiusError += 2 * (t - r) + 1;
		}
	}
};
template<typename T>
inline static void DrawSphere(vector<vector<vector<T>>> & grid, T value, int x0, int y0, int z0, int radius) {
	int r = radius, t = 0;
	int radiusError = 1 - r;

	while (r >= t) {
		// pass in base point (x0,y0,z0), this algorithm's t as the current plane,
		// this algorithm's r as the radius, and pass along radius error.
		DrawCircleXY(grid, value, x0, y0, z0 - t, r, radiusError);
		DrawCircleXY(grid, value, x0, y0, z0 + t, r, radiusError);
		DrawCircleYZ(grid, value, x0 - t, y0, z0, r, radiusError);
		DrawCircleYZ(grid, value, x0 + t, y0, z0, r, radiusError);
		DrawCircleXZ(grid, value, x0, y0 + t, z0, r, radiusError);
		DrawCircleXZ(grid, value, x0, y0 - t, z0, r, radiusError);
		++t;
		if (radiusError < 0) {
			radiusError += 2 * t + 1;
		} else {
			r--;
			radiusError += 2 * (t - r) + 1;
		}
	}
};
inline static void DrawBall(voxelGrid & grid, int x0, int y0, int z0, int radius) {
	int r = radius, t = 0;
	int radiusError = 1 - r;

	while (r >= t) {
		// pass in base point (x0,y0,z0), this algorithm's t as the current plane,
		// this algorithm's r as the radius, and pass along radius error.
		DrawDiskXY(grid, x0, y0, z0 - r, t, radiusError);
		DrawDiskXY(grid, x0, y0, z0 + r, t, radiusError);
		DrawDiskXY(grid, x0, y0, z0 - t, r, radiusError);
		DrawDiskXY(grid, x0, y0, z0 + t, r, radiusError);
		DrawCircleYZ(grid, x0 - t, y0, z0, r, radiusError);
		DrawCircleYZ(grid, x0 + t, y0, z0, r, radiusError);
		DrawCircleXZ(grid, x0, y0 + t, z0, r, radiusError);
		DrawCircleXZ(grid, x0, y0 - t, z0, r, radiusError);
		++t;
		if (radiusError < 0) {
			radiusError += 2 * t + 1;
		} else {
			r--;
			radiusError += 2 * (t - r) + 1;
		}
	}
};

#endif /* SPHERE_DRAWSPHERE_H_ */
