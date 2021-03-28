/*
 * normalizeGrid.h
 *
 *  Created on: Apr 19, 2015
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_NORMALIZEGRID_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_NORMALIZEGRID_H_

#include "../utils/numerical_limits.h"
#include <math.h>
#include <vector>

using namespace std;

/**
 * Simple grid normalization algorithm. The values of the input grid are mapped inside
 * the [0, 1] interval (min -> 0, max -> 1).
 *
 * @param grid	the grid containing the values to normalize
 */
template <typename T>
static inline void normalizeGrid(vector<vector<vector<T>>> & grid) {
	int length = grid.size();
	int width  = grid[0].size();
	int height = grid[0][0].size();

	T min = std::numeric_limits<T>::max();
	T max = -min;
	bool zero_grid = true;
	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if(grid[i][j][k] > 0) {
					zero_grid = false;
					if (min > grid[i][j][k])
					min = grid[i][j][k];
					if (max < grid[i][j][k])
					max = grid[i][j][k];
				}
			}
		}
	}
	T delta = max - min;
	if (zero_grid || fabs(delta) < 1.0e-10)
		return;

	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if(grid[i][j][k] > 0) {
					grid[i][j][k] -= min;
					grid[i][j][k] /= delta;
				}
			}
		}
	}
};
template <typename T>
static inline void normalizePositive(vector<vector<vector<T>>> & grid) {
	int length = grid.size();
	int width  = grid[0].size();
	int height = grid[0][0].size();

	T max = 0;

	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if (max < grid[i][j][k])
					max = grid[i][j][k];
			}
		}
	}
	if (max < 1.0e-20)
		return;

	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if (grid[i][j][k] > 0) {
					grid[i][j][k] /= max;
				}
			}
		}
	}
};
template <typename T>
static inline void normalizeNegative(vector<vector<vector<T>>> & grid) {
	int length = grid.size();
	int width  = grid[0].size();
	int height = grid[0][0].size();

	T min = 0;

	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if (min > grid[i][j][k])
					min = grid[i][j][k];
			}
		}
	}

	if (fabs(min) < 1.0e-20)
		return;

	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if (grid[i][j][k] < 0) {
					grid[i][j][k] /= fabs(min);
				}
			}
		}
	}
};

template <typename T>
static inline void standardiseGrid(vector<vector<vector<T>>> & grid) {
	int length = grid.size();
	int width  = grid[0].size();
	int height = grid[0][0].size();

	size_t N = 0;
	T mean = 0.0;
	T M2 = 0.0;

	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if(grid[i][j][k] != 0.0) {
					++N;
			        double delta = grid[i][j][k] - mean;
			        mean += delta / N;
			        M2 += delta * (grid[i][j][k] - mean);
				}
			}
		}
	}
	T stdev = sqrt(M2 / (N - 1));
	if (stdev < 1.0e-10)
		return;

	for (int i = 0; i < length; ++i) {
		for (int j = 0; j < width; ++j) {
			for (int k = 0; k < height; ++k) {
				if(grid[i][j][k] > 0) {
					grid[i][j][k] -= mean;
					grid[i][j][k] /= stdev;
				}
			}
		}
	}
};

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_NORMALIZEGRID_H_ */
