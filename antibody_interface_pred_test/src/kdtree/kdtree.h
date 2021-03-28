/*
 * kd-tree.h
 *
 *  Created on: 18/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_KDTREE_KDTREE_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_KDTREE_KDTREE_H_

#include "../SurfacePatch/CompactPatchDescriptor.h"
#include "nanoflann.hpp"
#include <vector>
#include "../Atom/atom.h"

using namespace std;
using namespace nanoflann;

typedef struct PointCloud_of_cpd {
	vector<CompactPatchDescriptor> const * points;
	/**
	 * Constructor
	 */
	inline PointCloud_of_cpd(vector<CompactPatchDescriptor> const * pts) {
		assert(pts != NULL);
		points = pts;
	};

	/**
	 * Must return the number of data points
	 */
	inline size_t kdtree_get_point_count() const {
		return points->size();
	};

	/**
	 * Returns the distance between the vector "p1[0:size-1]" and the
	 * data point with index "idx_p2" stored in the class:
	 */
	inline double kdtree_distance(const double *p1, const size_t idx_p2, size_t size) const {
		return (*points)[idx_p2].center.s_distance(point3D(p1[0], p1[1], p1[2]));
	};

	/**
	 * Returns the dim'th component of the idx'th point in the class:
	 * Since this is inlined and the "dim" argument is typically an immediate value,
	 * the "if/else's" are actually solved at compile time.
	 */
	inline double kdtree_get_pt(const size_t idx, int dim) const {
		if (dim == 0)
			return (*points)[idx].center.x;
		else if (dim == 1)
			return (*points)[idx].center.y;
		else
			return (*points)[idx].center.z;
	};
	/**
	 * Optional bounding-box computation: return false to default to a standard bbox computation loop.
	 * Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	 * Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	 */
	template<class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const {
		return false;
	};

} PointCloud_of_cpd;


typedef struct kdtree_cpd {
	typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud_of_cpd>, PointCloud_of_cpd, 3 /* dim */> kd_tree_cpd;

public:
	inline kdtree_cpd(vector<CompactPatchDescriptor> const * pts, size_t leaf_max_size = 2) {
		cloud = new PointCloud_of_cpd(pts);
		index = new kd_tree_cpd(3 /*dim*/, *cloud , KDTreeSingleIndexAdaptorParams(leaf_max_size));
		index->buildIndex();
	};
	virtual ~kdtree_cpd() {
		delete index;
		delete cloud;
	};
	/**
	 * Perform a search for the closest point.
	 * @param 	query			the query point
	 * @param	out_dist_sqr	the squared Euclidean distances of neighbors to the query point
	 * @param	eps				search for eps-approximate neighbors
	 *
	 * @return					the index of the search result
	 */
	inline size_t nnSearch(point3D const & query, double & out_dist_sqr, float eps = 0.0) {
		SearchParams params(32, eps, true);
		const double query_pt[3] = {query.x, query.y, query.z};
		size_t ret_index;

		KNNResultSet<double, size_t, size_t> resultSet(1);
		resultSet.init(&ret_index, &out_dist_sqr);
		index->findNeighbors(resultSet, &query_pt[0], params);
		return ret_index;
	};
	/**
	 * Perform a search for the N closest points.
	 * @param 	query			the query point
	 * @param	num_results		the maximum number of neighbors to return
	 * @param	ret_index		vector which will contain the indexes of the returned points
	 * @param	out_dist_sqr	vector which will contain the squared Euclidean distances of
	 * 							neighbors to the query point
	 * @param	eps				search for eps-approximate neighbors
	 *
	 */
	inline void knnSearch(point3D const & query, size_t num_results, vector<size_t> & ret_index,
			vector<double> & out_dist_sqr, float eps = 0.0) {
		SearchParams params(32, eps, true);

		const double query_pt[3] = {query.x, query.y, query.z};

		KNNResultSet<double, size_t, size_t> resultSet(num_results);
		resultSet.init(&ret_index[0], &out_dist_sqr[0]);
		index->findNeighbors(resultSet, &query_pt[0], params);

	};
	/**
	 * Find all the neighbors to a query point within a maximum radius
	 * @param 	query			the query point
	 * @param	search_radius	the radius of the search
	 * @param	ret_matches		The output is given as a vector of pairs, of which the first element
	 * 							is a point index and the second the corresponding distance.
	 *
	 * @param	eps				search for eps-approximate neighbors
	 * @param 	sorted			only for radius search, require neighbors sorted by distance
	 *
	 * @return					The number of points within the given radius (i.e. indices.size() or dists.size())
	 *
	 */
	inline size_t radiusSearch(point3D const & query, double search_radius,
			vector<std::pair<size_t, double> > & ret_matches, float eps = 0.0, bool sorted = true) {
		SearchParams params(32, eps, sorted);
		const double query_pt[3] = {query.x, query.y, query.z};
		return index->radiusSearch(&query_pt[0], search_radius, ret_matches, params);
	};
private:
	kd_tree_cpd * index;
	PointCloud_of_cpd * cloud;
} kdtree_cpd;


typedef struct PointCloud_of_atom {
	vector<atom> const * points;
	/**
	 * Constructor
	 */
	inline PointCloud_of_atom(vector<atom> const * pts) {
		assert(pts != NULL);
		points = pts;
	};

	inline PointCloud_of_atom(PointCloud_of_atom const & cloud) {
		assert(cloud.points != NULL);
		this->points = cloud.points;
	};

	/**
	 * Must return the number of data points
	 */
	inline size_t kdtree_get_point_count() const {
		return points->size();
	};

	/**
	 * Returns the distance between the vector "p1[0:size-1]" and the
	 * data point with index "idx_p2" stored in the class:
	 */
	inline float kdtree_distance(const float *p1, const size_t idx_p2, size_t size) const {
		float dx = (*points)[idx_p2].x - p1[0];
		float dy = (*points)[idx_p2].y - p1[1];
		float dz = (*points)[idx_p2].z - p1[2];
		return (dx*dx + dy*dy + dz*dz);
//		return sqrt(dx*dx + dy*dy + dz*dz);
	};

	/**
	 * Returns the dim'th component of the idx'th point in the class:
	 * Since this is inlined and the "dim" argument is typically an immediate value,
	 * the "if/else's" are actually solved at compile time.
	 */
	inline float kdtree_get_pt(const size_t idx, int dim) const {
		if (dim == 0)
			return (*points)[idx].x;
		else if (dim == 1)
			return (*points)[idx].y;
		else
			return (*points)[idx].z;
	};
	/**
	 * Optional bounding-box computation: return false to default to a standard bbox computation loop.
	 * Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	 * Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	 */
	template<class BBOX>
	bool kdtree_get_bbox(BBOX &bb) const {
		return false;
	};

} PointCloud_of_atom;


typedef struct kdtree_atom {
	typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<float, PointCloud_of_atom>, PointCloud_of_atom, 3 /* dim */> kd_tree_atom;

public:
	inline kdtree_atom(vector<atom> const * pts, size_t leaf_max_size = 2) {
		cloud = new PointCloud_of_atom(pts);
		leafMaxSize = leaf_max_size;
		index = new kd_tree_atom(3 /*dim*/, *cloud , KDTreeSingleIndexAdaptorParams(leaf_max_size));
		index->buildIndex();
	};
	inline kdtree_atom() : cloud(NULL), index(NULL), leafMaxSize(0) { };
	inline kdtree_atom(kdtree_atom const & tree) {
		this->cloud = new PointCloud_of_atom(*tree.cloud);
		this->leafMaxSize = tree.leafMaxSize;
		this->index = new kd_tree_atom(3 /*dim*/, *(this->cloud) , KDTreeSingleIndexAdaptorParams(this->leafMaxSize));
		this->index->buildIndex();

	};
	/**
	 * Copy assignment operator
	 */
	inline kdtree_atom & operator=(kdtree_atom const & tree) {
		// protect against invalid self-assignment
		if (this != &tree) {
			if (this->cloud != NULL)
				delete this->cloud;
			if (this->index != NULL)
				delete this->index;

			this->cloud = new PointCloud_of_atom(*tree.cloud);
			this->leafMaxSize = tree.leafMaxSize;
			this->index = new kd_tree_atom(3 /*dim*/, *(this->cloud) , KDTreeSingleIndexAdaptorParams(this->leafMaxSize));
			this->index->buildIndex();
		}
		// by convention, always return *this
		return *this;
	}
	virtual ~kdtree_atom() {
		delete index;
		delete cloud;
	};
	inline bool empty() const{
		return (cloud == NULL || cloud->points->empty());
	}
	/**
	 * Perform a search for the closest point.
	 * @param 	query			the query point
	 * @param	out_dist_sqr	the squared Euclidean distances of neighbors to the query point
	 * @param	eps				search for eps-approximate neighbors
	 *
	 * @return					the index of the search result
	 */
	inline size_t nnSearch(atom const & query, float & out_dist_sqr, float eps = 0.0) const {
		SearchParams params(32, eps, true);
		const float query_pt[3] = {query.x, query.y, query.z};
		size_t ret_index;

		KNNResultSet<float, size_t, size_t> resultSet(1);
		resultSet.init(&ret_index, &out_dist_sqr);
		index->findNeighbors(resultSet, &query_pt[0], params);
		return ret_index;
	};

	/**
	 * Perform a search for the closest point.
	 * @param 	query			the query point
	 * @param	out_dist_sqr	the squared Euclidean distances of neighbors to the query point
	 * @param	eps				search for eps-approximate neighbors
	 *
	 * @return					the index of the search result
	 */
	inline size_t nnSearch(point3D const & query, float & out_dist_sqr, float eps = 0.0) const {
		SearchParams params(32, eps, true);
		const float query_pt[3] = {(float)query.x, (float)query.y, (float)query.z};
		size_t ret_index;

		KNNResultSet<float, size_t, size_t> resultSet(1);
		resultSet.init(&ret_index, &out_dist_sqr);
		index->findNeighbors(resultSet, &query_pt[0], params);
		return ret_index;
	};
	/**
	 * Perform a search for the N closest points.
	 * @param 	query			the query point
	 * @param	num_results		the maximum number of neighbors to return
	 * @param	ret_index		vector which will contain the indexes of the returned points
	 * @param	out_dist_sqr	vector which will contain the squared Euclidean distances of
	 * 							neighbors to the query point
	 * @param	eps				search for eps-approximate neighbors
	 *
	 */
	inline void knnSearch(atom const & query, size_t num_results, vector<size_t> & ret_index,
			vector<float> & out_dist_sqr, float eps = 0.0) {
		SearchParams params(32, eps, true);

		const float query_pt[3] = {query.x, query.y, query.z};

		KNNResultSet<float, size_t, size_t> resultSet(num_results);
		resultSet.init(&ret_index[0], &out_dist_sqr[0]);
		index->findNeighbors(resultSet, &query_pt[0], params);

	};
	/**
	 * Find all the neighbors to a query point within a maximum radius
	 * @param 	query			the query point
	 * @param	search_radius	the radius of the search
	 * @param	ret_matches		The output is given as a vector of pairs, of which the first element
	 * 							is a point index and the second the corresponding distance.
	 *
	 * @param	eps				search for eps-approximate neighbors
	 * @param 	sorted			only for radius search, require neighbors sorted by distance
	 *
	 * @return					The number of points within the given radius (i.e. indices.size() or dists.size())
	 *
	 */
	inline size_t radiusSearch(atom const & query, float search_radius,
			vector<std::pair<size_t, float> > & ret_matches, float eps = 0.0, bool sorted = true) const {
		SearchParams params(32, eps, sorted);
		const float query_pt[3] = {query.x, query.y, query.z};
		return index->radiusSearch(&query_pt[0], search_radius, ret_matches, params);
	};
	inline size_t radiusSearch(point3D const & query, float search_radius,
			vector<std::pair<size_t, float> > & ret_matches, float eps = 0.0, bool sorted = true) const {
		SearchParams params(32, eps, sorted);
		const float query_pt[3] = {(float)query.x, (float)query.y, (float)query.z};
		return index->radiusSearch(&query_pt[0], search_radius, ret_matches, params);
	};
private:
	kd_tree_atom * index;
	PointCloud_of_atom * cloud;
	size_t leafMaxSize;
} kdtree_atom;
#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_KDTREE_KDTREE_H_ */
