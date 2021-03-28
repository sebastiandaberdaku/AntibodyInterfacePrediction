/*
 * SurfaceTriangle.h
 *
 *  Created on: 20/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_SURFACETRIANGLE_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_SURFACETRIANGLE_H_

#include "../SurfacePatch/CompactPatchDescriptor.h"
#include <unordered_map>

#define ANGLE_MAX (M_PI / 4)
#define DIST_MAX 4
#define DIST_MIN 1

using namespace std;

/**
 * There are six possible ways to position a triangle over another one, depending on
 * the chosen vertex association. Let's call the first triangle A1B1C1, where A1, B1
 * and C1 are the names of its vertices, and let's call the second triangle A2B2C2
 * where A2, B2 and C2 are its vertices. There are six possible combinations:
 * 1) {(A1, A2), (B1, B2), (C1, C2)} = ABC
 * 2) {(A1, A2), (B1, C2), (C1, B2)} = ACB
 * 3) {(A1, B2), (B1, A2), (C1, C2)} = BAC
 * 4) {(A1, B2), (B1, C2), (C1, A2)} = BCA
 * 5) {(A1, C2), (B1, A2), (C1, B2)} = CAB
 * 6) {(A1, C2), (B1, B2), (C1, A2)} = CBA
 */
typedef enum vertexPermutation { ABC, ACB, BAC, BCA, CAB, CBA } vertexPermutation;

/**
 * This struct defines a triangle formed by three distinct patch centers
 * on the same molecular surface. The patch normals must be coherent with
 * the triangle plane (within a certain tolerance), and the triangle
 * vertexes must be within a certain distance interval, i.e. AB >= DIST_MIN,
 * and AB <= DIST_MAX, where A and B are generic vertexes of a triangle.
 * In order to generate e docking pose, two surface triangles must match,
 * i.e. there must be high similarity between corresponding surface patches,
 * the patch centers must be within some tolerance distance, and the
 * patch normals of the first triangle must be coherent with those on the
 * second one.
 */
typedef struct SurfaceTriangle {
	CompactPatchDescriptor const *p1, *p2, *p3;
	SurfaceTriangle();
	SurfaceTriangle(CompactPatchDescriptor const & p1,
			CompactPatchDescriptor const & p2,
			CompactPatchDescriptor const & p3);
	SurfaceTriangle(CompactPatchDescriptor const * p1,
			CompactPatchDescriptor const * p2,
			CompactPatchDescriptor const * p3);
	SurfaceTriangle(SurfaceTriangle const & tr);
	/**
	 * Copy assignment operator
	 */
	SurfaceTriangle & operator=(SurfaceTriangle const & tr);
	double getPerimeter();
	double getSArea();
	double getArea();
	friend ostream & operator<<(ostream & os, SurfaceTriangle const & t);

	static bool compatibleSurfaceTriangle(CompactPatchDescriptor const & p1,
			CompactPatchDescriptor const & p2, CompactPatchDescriptor const & p3);
	static void generateSurfaceTriangles(vector<CompactPatchDescriptor> const & patchList,
			vector<SurfaceTriangle> & surfTriangles);
	static void generateSurfaceTriangles(vector<CompactPatchDescriptor> const & patchList,
			int k, vector<SurfaceTriangle> & surfTriangles);
	static void generateSurfaceTriangles(vector<CompactPatchDescriptor> const & patchList,
			float radius, vector<SurfaceTriangle> & surfTriangles);
	static void createTriangleMultimap(vector<SurfaceTriangle> const & surfTriangleList,
			unordered_multimap<CompactPatchDescriptor const *, SurfaceTriangle const *> & surfTriangleMap);
	SurfaceTriangle trPermutation1(CompactPatchDescriptor const * cpd) const;
	SurfaceTriangle trPermutation2(CompactPatchDescriptor const * cpd) const;

	double getSimilarityScore(SurfaceTriangle const & t, vertexPermutation alignment = ABC) const;

private:
	/**
	 * Simple method to calculate Pearson's correlation between two
	 * three-dimensional arrays.
	 * @param x		first data array
	 * @param y		second data array
	 * @return		the Pearson's coefficient, a number in [-1, 1]
	 */
	static inline double correlation(double x[3], double y[3]) {
		double Ex = (x[0] + x[1] + x[2]) / 3;
		double Ey = (y[0] + y[1] + y[2]) / 3;

		double Sxx = 0, Syy = 0, Sxy = 0;

		for (int ii = 0; ii < 3; ++ii) {
			double dx = x[ii] - Ex;
			double dy = y[ii] - Ey;
			Sxx += dx * dx;
			Syy += dy * dy;
			Sxy += dx * dy;
		}
		/*
		 * TINY_VALUE could be something like 1e-20 and is used to "compensate"
		 * perfect correlation case (and avoid special verification).
		 */
		double TINY_VALUE = 1e-20;
		return Sxy / (sqrt(Sxx * Syy) + TINY_VALUE);
	};

} SurfaceTriangle;

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_SURFACETRIANGLE_H_ */
