/*
 * ZernikeDescriptor.cpp
 *
 *  Created on: Jul 8, 2014
 *      Author: sebastian
 */

#include "../utils/makeDirectory.h"
#include "ZernikeDescriptor.h"



ZernikeDescriptor::ZernikeDescriptor(array3D const & voxels, point3D const & center,
		double scale, int order, bool normalize) :
	order(order), voxels(&voxels), COG(center), scale(scale), length(voxels.size()),
	width(voxels[0].size()), height(voxels[0][0].size()) {
	assert(order > 0);

	computeGeometricMoments();

	if (fabs(gm->getMoment(0, 0, 0)) > 0) {
		computeZernikeMoments();
		computeInvariants();
		if (normalize)
			normalizeInvariants();
	} else {
//		cout << "WARNING: Computing zero invariants.\n";
		zm = NULL;
		gm = NULL;
		computeZeroInvariants();
	}
}

ZernikeDescriptor::ZernikeDescriptor() : voxels(NULL), gm(NULL), zm(NULL), length(0), width(0), height(0), order(0), scale(0) {}

ZernikeDescriptor::~ZernikeDescriptor() {
	if (zm != NULL)
		delete zm;
	if (gm != NULL)
		delete gm;
}
/**
 * Computes the Zernike moments.
 */
void ZernikeDescriptor::computeZernikeMoments() {
	// Zernike moments
	zm = new ZernikeMoments(order, gm);
}
/**
 * Computes the geometric moments.
 */
void ZernikeDescriptor::computeGeometricMoments() {
	// Geometric moments
	gm = new GeometricMoments(*voxels, scale, COG, order);
}
/**
 * Computes the geometric moments first, then the Zernike moments.
 */
void ZernikeDescriptor::computeZeroMoments() {
	// Geometric moments
	gm = new GeometricMoments();
	// Zernike moments
	zm = new ZernikeMoments();
}
/**
 * Reconstructs the original object from the 3D Zernike moments.
 * @param grid		result grid as 3D complex stl vector
 * @param minN		min value for n freq index
 * @param maxN 		max value for n freq index
 * @param minL		min value for l freq index
 * @param maxL		max value for l freq index
 */
void ZernikeDescriptor::reconstruct(complexArray3D& grid, int minN, int maxN, int minL, int maxL) {
	zm->reconstruct(grid, COG, scale);
}

/**
 * Computes the Zernike moment based invariants, i.e. the norms of vectors with
 * components of Z_nl^m with m being the running index.
 */

void ZernikeDescriptor::computeInvariants() {
	invariants.clear();
	for (int n = 0; n <= order; ++n) {
		for (int l = n % 2; l <= n; l += 2) {
			double sum = 0;
			for (int m = -l; m <= l; ++m) {
				/**
				 * The norm value of a complex number is its squared magnitude,
				 * defined as the addition of the square of both its real and
				 * its imaginary part (without the imaginary unit). This is
				 * the square of abs(x).
				 */
				sum += norm(zm->getMoment(n, l, m));
			}
			invariants.push_back(sqrt(sum));
		}
	}
//	double zmi = invariants[0];
//	if (zmi > 0)
//		for (auto & inv : invariants)
//			inv /= zmi;
}
/**
 * Computes the Zernike zero invariants.
 */

void ZernikeDescriptor::computeZeroInvariants() {
	invariants.clear();
	for (int n = 0; n <= order; ++n) {
		for (int l = n % 2; l <= n; l += 2) {
			invariants.push_back(0);
		}
	}
}
/**
 * Saves the computed invariants into a binary file
 * @param fName		name of the output file
 */
void ZernikeDescriptor::saveInvariants(string const & fName) {
	ofstream file_stream;
	makeDirectory("./output");
	file_stream.open(("./output/" + fName + ".inv").c_str());
	file_stream.precision(6);
	file_stream.setf(ios::scientific, ios::floatfield);
	int dim = invariants.size();
	file_stream << "Number of Invariants: " << dim << "\n";

	for (int i = 0; i < dim; ++i) {
		file_stream << invariants[i] << "\t";
	}
	file_stream.close();
}
