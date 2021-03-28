/*
 * ZernikeDescriptor.h
 *
 *  Created on: Jul 8, 2014
 *      Author: sebastian
 */

#ifndef ZERNIKEDESCRIPTOR_H_
#define ZERNIKEDESCRIPTOR_H_

#include "../Geometry/point3D.h"
#include "../Geometry/voxelGrid.h"
#include "GeometricMoments.h"
#include "ZernikeMoments.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

/**
 * This class serves as a wrapper around the geometric and
 * Zernike moments. It provides also the implementation of invariant Zernike
 * descriptors, means of reconstruction of orig. function, etc.
 */
typedef class ZernikeDescriptor {
public:
	/**
	 * Constructor of the class.
	 * @param voxels 	pointer to the voxel grid containing the 3D function
	 * @param order		maximal order of the Zernike moments (N in paper)
	 */
	ZernikeDescriptor(array3D const & voxels, point3D const & center, double scale, int order, bool normalize = false);
	/**
	 * Default constructor.
	 */
	ZernikeDescriptor();
	virtual ~ZernikeDescriptor();
	void reconstruct(complexArray3D& grid, int minN = 0, int maxN = 100,
			int minL = 0, int maxL = 100);
	void saveInvariants(string const & fName);
    array1D invariants;			//vector of invariants under SO(3)
    /**
     * Extraction operator
     * The overload of operator<< that takes a std::ostream& as the left hand argument.
     * Since this operator takes the user-defined type as the right argument (b in a@b),
	 * it must be implemented as non-members.
     */
	friend std::ostream & operator<<(std::ostream& os, ZernikeDescriptor const & zd) {
//		os.precision(std::numeric_limits<double>::digits10 + 2);
//		os.setf(ios::scientific, ios::floatfield);
		int size = zd.invariants.size();
		for(int i = 0; i < size; ++i) {
			os << zd.invariants[i] << ", ";
		}
		os << "\n";
		return os;
	};
	void printGeometricMoments() {
		for (int r = 0; r <= order; ++r) {
			for (int s = 0; r + s <= order; ++s) {
				for (int t = 0; r + s + t <= order; ++t) {
					cout << "GM[" << r << ", " << s << ", " << t << "] = "
							<< gm->getMoment(r, s, t) << "; ";
				}
			}
		}
		cout << endl;
	};
	void printZernikeMoments() {
		for (int n = 0; n <= order; ++n) {
			for (int l = n % 2; l <= n; l += 2) {
				for (int m = 0; m <= l; ++m) {
					cout << "ZM[" << n << ", " << l << ", " << m << "] = "
							<< zm->getMoment(n, l, m) << "; ";
				}
			}
		}
		cout << endl;
	};
private:
	/**
	 * The 3D Zernike descriptor is mathematically most noteworthy rotationally invariant.
	 * This is one of the largest advantages of 3D Zernike descriptor. However, in practice
	 * the descriptor of rotated protein structures are not exactly identical. This error
	 * is caused possibly when the protein surface shape is discretized into voxels. We
	 * found that in computing all the distance measures, that is, Euclidean, Manhattan,
	 * and the correlation coefficient-based, normalizing each number in a 3D Zernike
	 * descriptor by the sum of the 121 numbers of the descriptor reduces the error the best
	 * among tested methods.
	 */
	inline void normalizeInvariants() {
		double norm = 0;
		for (auto & i : invariants)
			norm += i * i;

		if (norm > 0)
			for (auto & i : invariants)
				i /= sqrt(norm);
	};
	// ---- private helper functions ----
	void computeGeometricMoments();
	void computeZernikeMoments();
	void computeZeroMoments();
	void computeInvariants();
	void computeZeroInvariants();
	/**
	 * Determines the scaling factor in order to scale and translate our 3D surface
	 * inside the unit ball.
	 * @param COG		the center of gravity of the 3D function
	 * @return			the scaling factor
	 */
	inline double scaling_BuondingSphere(point3D const & COG) {
		double max = 0; /* maximum distance from the COG */
		double s_distance = 0;
		for (int i = 0; i < length; ++i) {
			for (int j = 0; j < width; ++j) {
				for (int k = 0; k < height; ++k) {
					if (fabs((*voxels)[i][j][k]) > 0) {
						s_distance = point3D(i, j, k).s_distance(COG);
						if (max < s_distance)
							max = s_distance;
					}
				}
			}
		}
		if (max < 1.0e-10)
			return 1.0;

		return 1 / (sqrt(max) + sqrt(3.0));
	};

	inline double scaling_RadiusVar(point3D const & COG) {
		double sum = 0; /* maximum distance from the COG */
		int nVox = 0;
		for (int i = 0; i < length; ++i) {
			for (int j = 0; j < width; ++j) {
				for (int k = 0; k < height; ++k) {
					if (fabs((*voxels)[i][j][k]) > 0) {
						sum += point3D(i, j, k).s_distance(COG);
						++nVox;
					}
				}
			}
		}
		if (fabs(sum) < 1.0e-10)
			return 1.0;

		return 1 / (2 * sqrt(sum/nVox));
	}

private:
    // ---- member variables ----
    int     order;              	// maximal order of the moments to be computed (max{n})
    int     length, width, height;  // dimensions of the voxel grid
    double 	scale;
    array3D const * voxels;			// 3D array containing the voxels
    point3D COG;             		// center of gravity
public:
    ZernikeMoments *     zm;
    GeometricMoments *   gm;

} ZernikeDescriptor;
#endif /* ZERNIKEDESCRIPTOR_H_ */
