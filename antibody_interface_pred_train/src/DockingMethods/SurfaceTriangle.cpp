///*
// * SurfaceTriangle.cpp
// *
// *  Created on: 25/feb/2015
// *      Author: sebastian
// */
//
//#include "../kdtree/kdtree.h"
//#include "SurfaceTriangle.h"
//
//SurfaceTriangle::SurfaceTriangle() :
//		p1(NULL), p2(NULL), p3(NULL) {
//}
//
//SurfaceTriangle::SurfaceTriangle(CompactPatchDescriptor const & p1,
//		CompactPatchDescriptor const & p2, CompactPatchDescriptor const & p3) :
//		p1(&p1), p2(&p2), p3(&p3) {
//}
//
//SurfaceTriangle::SurfaceTriangle(CompactPatchDescriptor const * p1,
//		CompactPatchDescriptor const * p2, CompactPatchDescriptor const * p3) :
//		p1(p1), p2(p2), p3(p3) {
//}
//
//SurfaceTriangle::SurfaceTriangle(SurfaceTriangle const & tr) :
//		p1(tr.p1), p2(tr.p2), p3(tr.p3) {
//}
//SurfaceTriangle & SurfaceTriangle::operator=(SurfaceTriangle const & tr) {
//	if (this != &tr) {
//		this->p1 = tr.p1;
//		this->p2 = tr.p2;
//		this->p3 = tr.p3;
//	}
//	return *this;
//}
//
//double SurfaceTriangle::getPerimeter() {
//	double a = p1->center.distance(p2->center);
//	double b = p1->center.distance(p3->center);
//	double c = p2->center.distance(p3->center);
//	return a + b + c;
//}
///**
// * The shape of the triangle is determined by the lengths
// * of the sides. Therefore the area can also be derived
// * from the lengths of the sides. By Heron's formula:
// *
// * $A = \sqrt{s(s-a)(s-b)(s-c)}$
// *
// * where $s = \frac{a+b+c}{2}$ is the semiperimeter, or
// * half of the triangle's perimeter.
// */
//double SurfaceTriangle::getSArea() {
//	double a = p1->center.distance(p2->center);
//	double b = p1->center.distance(p3->center);
//	double c = p2->center.distance(p3->center);
//	double s = (a + b + c) / 2;
//	return s * (s - a) * (s - b) * (s - c);
//}
//
//double SurfaceTriangle::getArea() {
//	return sqrt(this->getSArea());
//}
//
//ostream & operator<<(ostream & os, SurfaceTriangle const & t) {
//	os << "Triangle: [" << t.p1->center << ", " << t.p2->center << ", "
//			<< t.p3->center << "]\n";
//	os << "AB: " << t.p1->center.distance(t.p2->center);
//	os << " BC: " << t.p2->center.distance(t.p3->center);
//	os << " AC: " << t.p1->center.distance(t.p3->center) << "\n";
//
////		os << "Descriptor 1: " << t.p1;
////		os << "Descriptor 2: " << t.p2;
////		os << "Descriptor 3: " << t.p3 << "\n";
//	return os;
//}
//
///**
// * This method checks if the three input surface patches can be used to form
// * a SurfaceTriangle, by checking the distance and angle constraints.
// * Returns true if a SurfaceTriangle can be created with the three input CPD.
// */
//bool SurfaceTriangle::compatibleSurfaceTriangle(CompactPatchDescriptor const & p1,
//		CompactPatchDescriptor const & p2, CompactPatchDescriptor const & p3) {
//
//	if (p1.ID >= p2.ID || p2.ID >= p3.ID)
//		return false;
//
//	point3D A = p1.center;
//	point3D B = p2.center;
//	point3D C = p3.center;
//
//	double bc = B.distance(C);
//	double ac = A.distance(C);
//	double ab = A.distance(B);
//
//	if (bc > DIST_MAX || bc < DIST_MIN
//			|| ac > DIST_MAX || ac < DIST_MIN
//			|| ab > DIST_MAX || ab < DIST_MIN
//			|| bc + ac <= ab || bc + ab <= ac || ac + ab <= bc)
//		return false;
//
//	point3D nA = p1.patchNormal;
//	point3D nB = p2.patchNormal;
//	point3D nC = p3.patchNormal;
//
//	if (nA.angle(nB) > ANGLE_MAX || nA.angle(nC) > ANGLE_MAX
//			|| nB.angle(nC) > ANGLE_MAX)
//		return false;
//
//	point3D AB = B - A;
//	point3D AC = C - A;
//	point3D nABC = AB * AC; /**< The normal vector to the ABC triangle plane. */
//
//	double a = nA.angle(nABC);
//	double b = nB.angle(nABC);
//	double c = nC.angle(nABC);
//
//	if ((a > ANGLE_MAX && a < M_PI - ANGLE_MAX)
//			|| (b > ANGLE_MAX && b < M_PI - ANGLE_MAX)
//			|| (c > ANGLE_MAX && c < M_PI - ANGLE_MAX))
//		return false;
//
//	return true;
//}
//
///**
// * This method generates SurfaceTriangles from a given list of CompactPatchDescriptors.
// * In order to ensure that the SurfaceTriangles are generated only once, triangles are
// * generated as ordered triples of distinct CompactPatchDescriptors.
// *
// * @param patchList		input list of CompactPatchDescriptors from the current
// * 						molecular surface
// * @param surfTriangles	output list of generated SurfaceTriangles
// */
//void SurfaceTriangle::generateSurfaceTriangles(vector<CompactPatchDescriptor> const & patchList,
//		vector<SurfaceTriangle> & surfTriangles) {
//	size_t n = patchList.size();
//	assert(n>=3);
//	for (size_t i = 0; i < n; ++i) {
//		for (size_t j = i + 1; j < n; ++j) {
//			for (size_t k = j + 1; k < n; ++k) {
//				if (compatibleSurfaceTriangle(patchList[i], patchList[j], patchList[k])) {
//					surfTriangles.push_back(SurfaceTriangle(patchList[i], patchList[j], patchList[k]));
//				}
//			}
//		}
//	}
//}
//
///**
// * This method generates SurfaceTriangles from a given list of CompactPatchDescriptors.
// * A kd-tree containing the centers of each surface patch is created. For each surface patch
// * center, the k-nearest neighbors are extracted from the kd-tree and triangles are formed
// * using those neighbors.
// *
// * @param patchList		input list of CompactPatchDescriptors from the current
// * 						molecular surface
// * @param k				number of neighbors
// * @param surfTriangles	output list of generated SurfaceTriangles
// */
//void SurfaceTriangle::generateSurfaceTriangles(vector<CompactPatchDescriptor> const & patchList,
//		int k, vector<SurfaceTriangle> & surfTriangles) {
//	size_t n = patchList.size();
//	assert(n >= k + 1);
//
//	vector<size_t> ret_idx;
//	vector<double> ret_sqr_dist;
//	ret_idx.resize(k);
//	ret_sqr_dist.resize(k);
//
//	kdtree_cpd cpdTree(&patchList);
//	for (size_t i = 0; i < n; ++i) {
//		CompactPatchDescriptor const * const cpd1 = &patchList[i];
//
//		cpdTree.knnSearch(cpd1->center, k, ret_idx, ret_sqr_dist);
//		assert(doubleCompare(ret_sqr_dist[0], 0));
//		for (size_t j = 1; j < k; ++j) {
//			CompactPatchDescriptor const * const cpd2 = &patchList[ret_idx[j]];
//			for (size_t l = 1; l < k; ++l) {
//				CompactPatchDescriptor const * const cpd3 = &patchList[ret_idx[l]];
//				if (compatibleSurfaceTriangle(*cpd1, *cpd2, *cpd3)) {
//					surfTriangles.push_back(SurfaceTriangle(cpd1, cpd2, cpd3));
//				}
//			}
//		}
//	}
//}
//
//
///**
// * This method generates SurfaceTriangles from a given list of CompactPatchDescriptors.
// * A kd-tree containing the centers of each surface patch is created. For each surface patch
// * center, the neighbors within a certain distance are extracted from the kd-tree and triangles
// * are formed using those neighbors.
// *
// * @param patchList		input list of CompactPatchDescriptors from the current
// * 						molecular surface
// * @param radius		radius of the search
// * @param surfTriangles	output list of generated SurfaceTriangles
// */
//void SurfaceTriangle::generateSurfaceTriangles(vector<CompactPatchDescriptor> const & patchList,
//		float radius, vector<SurfaceTriangle> & surfTriangles) {
//	size_t n = patchList.size();
//	assert(n >= 3);
//
//	vector<std::pair<size_t, double> > ret_matches;
//
//	kdtree_cpd cpdTree(&patchList);
//#pragma omp parallel for private(ret_matches)
//	for (size_t i = 0; i < n; ++i) {
//		CompactPatchDescriptor const * const cpd1 = &patchList[i];
//
//		size_t k = cpdTree.radiusSearch(cpd1->center, radius, ret_matches);
//		assert(doubleCompare(ret_matches[0].second, 0));
//		for (size_t j = 1; j < k; ++j) {
//			CompactPatchDescriptor const * const cpd2 = &patchList[ret_matches[j].first];
//			if (cpd2->ID <= cpd1->ID)
//				continue;
//			for (size_t l = 1; l < k; ++l) {
//				CompactPatchDescriptor const * const cpd3 = &patchList[ret_matches[l].first];
//				if (cpd3->ID <= cpd2->ID)
//					continue;
//				if (compatibleSurfaceTriangle(*cpd1, *cpd2, *cpd3)) {
//#pragma omp critical
//					surfTriangles.push_back(SurfaceTriangle(cpd1, cpd2, cpd3));
//				}
//			}
//		}
//	}
//}
//
///**
// * Simple method that generates a SurfaceTriangle multimap. Each SurfaceTriangle is a
// * triple (a, b, c) of distinct CompactPatchDescriptors (pointers to CPDs). Given the
// * ID of a CompactPatchDescriptor, we want to know all the SurfaceTriangles in which
// * it appears. Thus, we use the CPD's ID as a key for the multimap, and the corresponding
// * SurfaceTriangle's IDs as values.
// *
// * @param surfTriangleList		input SurfaceTriangles list
// * @param surfTriangleMap		output multimap
// */
//void SurfaceTriangle::createTriangleMultimap(vector<SurfaceTriangle> const & surfTriangleList,
//		unordered_multimap<CompactPatchDescriptor const *, SurfaceTriangle const *> & surfTriangleMap) {
//	size_t n = surfTriangleList.size();
//	if (!surfTriangleMap.size())
//		surfTriangleMap.clear();
//	surfTriangleMap.reserve(3 * n);
//
//	using map_pair = pair<CompactPatchDescriptor const *, SurfaceTriangle const *>;
//
//	/**
//	 * The given triangle at index i appears is inserted as value for each of the three
//	 * IDs of its CPDs.
//	 */
//	for(size_t i = 0; i < n; ++i) {
//		surfTriangleMap.insert(map_pair(surfTriangleList[i].p1, &surfTriangleList[i]));
//		surfTriangleMap.insert(map_pair(surfTriangleList[i].p2, &surfTriangleList[i]));
//		surfTriangleMap.insert(map_pair(surfTriangleList[i].p3, &surfTriangleList[i]));
//	}
//}
//
///**
// * Permutes the input SurfaceTriangles in  such a way that the first CPD
// * of the permuted triangle equals the specified input CPD. There are two
// * possible permutations, this method only calculates one. The other permutation
// * is given by the next method for better performance.
// */
//SurfaceTriangle SurfaceTriangle::trPermutation1(CompactPatchDescriptor const * cpd) const {
//	if (cpd == this->p1) {
//		return SurfaceTriangle(this->p1, this->p2, this->p3);
//	} else if (cpd == this->p2) {
//		return SurfaceTriangle(this->p2, this->p1, this->p3);
//	} else if (cpd == this->p3) {
//		return SurfaceTriangle(this->p3, this->p2, this->p1);
//	} else {
//		throw invalid_argument("DockingPose::trPermutation1() - Unknown CompactPatchDescriptor ID!");
//		return SurfaceTriangle();
//	}
//}
///**
// * Permutes the input SurfaceTriangles in  such a way that the first CPD's ID
// * of the permuted triangle equals the specified input ID. There are two
// * possible permutations, this method only calculates one. The other permutation
// * is given by the previous method for better performance.
// */
//SurfaceTriangle SurfaceTriangle::trPermutation2(CompactPatchDescriptor const * cpd) const{
//	if (cpd == this->p1) {
//		return SurfaceTriangle(this->p1, this->p3, this->p2);
//	} else if (cpd == this->p2) {
//		return SurfaceTriangle(this->p2, this->p3, this->p1);
//	} else if (cpd == this->p3) {
//		return SurfaceTriangle(this->p3, this->p1, this->p2);
//	} else {
//		throw invalid_argument("DockingPose::trPermutation2() - Unknown CompactPatchDescriptor ID!");
//		return SurfaceTriangle();
//	}
//}
//
///**
// * This method calculates the similarity score between the current SurfaceTriangle
// * and the input one.
// * See the following code for details.
// *
// * @param t				input SurfaceTriangle
// * @param alignment		the vertexPermutation for the current match
// * @return				a measure of similarity between the SurfaceTriangles (high values indicate
// * 						high similarity)
// */
//double SurfaceTriangle::getSimilarityScore(SurfaceTriangle const & t, vertexPermutation alignment) const {
//	/** vertices of the first triangle */
//	point3D A1(this->p1->center), B1(this->p2->center), C1(this->p3->center);
//	/** normals on the vertices of the first triangle */
//	point3D nA1(this->p1->patchNormal), nB1(this->p2->patchNormal), nC1(this->p3->patchNormal);
//	/** vertices of the second triangle */
//	point3D A2, B2, C2;
//	/** normals on the vertices of the second triangle */
//	point3D nA2, nB2, nC2;
//	switch (alignment) {
//	case ABC:
//		A2 = t.p1->center; B2 = t.p2->center; C2 = t.p3->center;
//		nA2 = t.p1->patchNormal; nB2 = t.p2->patchNormal; nC2 = t.p3->patchNormal;
//		break;
//	case ACB:
//		A2 = t.p1->center; B2 = t.p3->center; C2 = t.p2->center;
//		nA2 = t.p1->patchNormal; nB2 = t.p3->patchNormal; nC2 = t.p2->patchNormal;
//		break;
//	case BAC:
//		A2 = t.p2->center; B2 = t.p1->center; C2 = t.p3->center;
//		nA2 = t.p2->patchNormal; nB2 = t.p1->patchNormal; nC2 = t.p3->patchNormal;
//		break;
//	case BCA:
//		A2 = t.p2->center; B2 = t.p3->center; C2 = t.p1->center;
//		nA2 = t.p2->patchNormal; nB2 = t.p3->patchNormal; nC2 = t.p1->patchNormal;
//		break;
//	case CAB:
//		A2 = t.p3->center; B2 = t.p1->center; C2 = t.p2->center;
//		nA2 = t.p3->patchNormal; nB2 = t.p1->patchNormal; nC2 = t.p2->patchNormal;
//		break;
//	case CBA:
//		A2 = t.p3->center; B2 = t.p2->center; C2 = t.p1->center;
//		nA2 = t.p3->patchNormal; nB2 = t.p2->patchNormal; nC2 = t.p1->patchNormal;
//		break;
//	default:
//		throw invalid_argument("DockingPose::matchSurfaceTriangles() - Unknown vertex permutation!");
//		break;
//	}
//
//	/* edge length */
//	double ab1 = A1.distance(B1), ac1 = A1.distance(C1), bc1 = B1.distance(C1);
//	double ab2 = A2.distance(B2), ac2 = A2.distance(C2), bc2 = B2.distance(C2);
//
//	/* angles between patch normals */
//	double nAB1 = nA1.angle(nB1), nAC1 = nA1.angle(nC1), nBC1 = nB1.angle(nC1);
//	double nAB2 = nA2.angle(nB2), nAC2 = nA2.angle(nC2), nBC2 = nB2.angle(nC2);
//
//	/**
//	 * Array containing the lengths of the edges for the two triangles.
//	 * Perfect correlation between these two arrays implies that the
//	 * triangles are equal.
//	 */
//	double edges1[3] = {ab1, ac1, bc1};
//	double edges2[3] = {ab2, ac2, bc2};
//
//	/**
//	 * Array containing the angles between patch normals.
//	 * Perfect correlation between these two arrays implies that the
//	 * normals are perfectly aligned.
//	 */
//	double angles1[3] = {nAB1, nAC1, nBC1};
//	double angles2[3] = {nAB2, nAC2, nBC2};
//
//	double triangle_shape_correlation = fabs(correlation(edges1, edges2));
//	double normals_alignment_correlation = fabs(correlation(angles1, angles2));
//	/*
//	 * High correlation: 0.5 to 1.0 or -0.5 to -1.0
//	 * Medium correlation: 0.3 to 0.5 or -0.3 to -0.5
//	 * Low correlation: 0.1 to 0.3 or -0.1 to -0.3
//	 */
//	if (triangle_shape_correlation < 0.3 || normals_alignment_correlation < 0.3)
//		return min_double;
//
//	//TODO da sistemare!!!
//
//	double descriptor_score = 0;
////	if (this->p1->type * t.p1->type <= 0)
////		descriptor_score += 1.0 + this->p1->inv_correlation(*t.p1);
////	if (this->p2->type * t.p2->type <= 0)
////		descriptor_score += 1.0 + this->p2->inv_correlation(*t.p2);
////	if (this->p3->type * t.p3->type <= 0)
////		descriptor_score += 1.0 + this->p3->inv_correlation(*t.p3);
//
//	return descriptor_score;
//}
