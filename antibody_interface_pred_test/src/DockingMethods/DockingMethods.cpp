/*
 * DockingMethods.cpp
 *
 *  Created on: Feb 26, 2015
 *      Author: sebastian
 */


#include "DockingMethods.h"
#include <assert.h>
/**
 * This method evaluates the current pose and relative score and decides whether to keep the pose or not.
 * A maximum number of allowed docking poses is provided. If the maximum number is reached, the method
 * evaluates if the current pose's score is greater than the worst score present in the list. If the
 * current score is smaller or equal than the worst score, the current pose is discarded, otherwise the
 * current pose is inserted in the poses list and the worst pose is discarded. The number of poses must
 * not exceed max_num_poses.
 *
 * @param currentPose	current pose to be evaluated
 * @param currentScore	current pose's score
 * @param max_num_poses	maximum number of poses
 * @param dockingPoses_map		map data-structure used to enforce the uniqueness of the poses in the list
 * @param orderedDockingPoses	ordered set data-structure used to mantain the poses ordered by
 * 								their scores
 */
//void DockingMethods::addPose(DockingPose const & currentPose, double currentScore, size_t max_num_poses,
//		unordered_map<DockingPose, double> & dockingPoses_map, set<scoredDockingPose, scoredDockingPose_compare> & orderedDockingPoses) {
//	/*
//	 * If the number of poses currently in the list
//	 * has reached its maximum.
//	 */
//	if (dockingPoses_map.size() == max_num_poses) {
//		/*
//		 * If the current pose has a better score than the worst one
//		 * currently in the list, find the correct position where
//		 * to insert it in order to keep the list ordered.
//		 */
//		if (currentScore > orderedDockingPoses.begin()->second) {
//			auto oldPose = dockingPoses_map.find(currentPose);
//			if (oldPose == dockingPoses_map.end()) {
//				// no duplicates
//				dockingPoses_map.insert(make_pair(currentPose, currentScore)); // insert currentPose
//				orderedDockingPoses.insert(make_pair(currentPose, currentScore)); // insert current pose in the ordered set
//				// remove the element with the least score
//				auto worst_element = orderedDockingPoses.begin();
//				dockingPoses_map.erase(worst_element->first);
//				orderedDockingPoses.erase(worst_element);
//			} else {
//				// existing duplicate found
//				if (currentScore <= oldPose->second) // if currentScore <= oldScore do nothing
//					return;
//				// erase old value
//				orderedDockingPoses.erase(make_pair(oldPose->first, oldPose->second));
//				dockingPoses_map.erase(oldPose);
//				// insert new value
//				dockingPoses_map.insert(make_pair(currentPose, currentScore)); // insert currentPose
//				orderedDockingPoses.insert(make_pair(currentPose, currentScore)); // insert current pose in the ordered set
//			}
//		}
//	}
//	/*
//	 * If there are less than max_num_poses poses push the current poses in
//	 * the list if not already present.
//	 */
//	else {
//		auto oldPose = dockingPoses_map.find(currentPose);
//		if (oldPose == dockingPoses_map.end()) {
//			// no duplicates
//			dockingPoses_map.insert(make_pair(currentPose, currentScore)); // insert currentPose
//			orderedDockingPoses.insert(make_pair(currentPose, currentScore)); // insert current pose in the ordered set
//		} else { // existing duplicate found
//			if (currentScore <= oldPose->second) // if currentScore <= oldScore do nothing
//				return;
//			// erase old value
//			orderedDockingPoses.erase(make_pair(oldPose->first, oldPose->second));
//			dockingPoses_map.erase(oldPose);
//			// insert new value
//			dockingPoses_map.insert(make_pair(currentPose, currentScore)); // insert currentPose
//			orderedDockingPoses.insert(make_pair(currentPose, currentScore)); // insert current pose in the ordered set
//		}
//	}
//};

/**
 * This method calculates and outputs the max_num_poses best DockingPoses.
 * The output vector will contain at most max_num_poses elements, ordered
 * in descending order by the quality of the pose, given by the score entry
 * in DockingPose.
 * @param bestPairs			input vector containing the best matching CPDpairs
 * @param surfTriangleMap_r	input multimap of the SurfaceTriangles for the
 * 							first molecular surface
 * @param surfTriangleMap_l	input multimap of the SurfaceTriangles for the
 * 							second molecular surface
 * @param max_num_poses		maximum number of DockingPoses to produce
 * @param dockingPoses		output vector containing the num_poses best
 * 							DockingPoses calculated by the method
 * @param dockingScores		list of poses ordered by score
 */
//void DockingMethods::calculateBestComplementarityPoses(vector<CPDpair> const & bestPairs,
//			unordered_multimap<CompactPatchDescriptor const *, SurfaceTriangle const *> & surfTriangleMap_r,
//			unordered_multimap<CompactPatchDescriptor const *, SurfaceTriangle const *> & surfTriangleMap_l,
//			size_t max_num_poses, vector<DockingPose> & dockingPoses, orderedPosesList & dockingScores) {
//	if (!dockingPoses.empty())
//		dockingPoses.clear();
//	if (!dockingScores.scores_list.empty())
//		dockingScores.scores_list.clear();
//	dockingPoses.reserve(max_num_poses); /**< reserve the space for max_num_poses entries */
//	dockingScores.reserve(max_num_poses);
//
//	unordered_map<DockingPose, double> dockingPoses_map;
//	dockingPoses_map.reserve(max_num_poses);
//	set<scoredDockingPose, scoredDockingPose_compare> orderedDockingPoses;
//
//	using m_it = unordered_multimap<CompactPatchDescriptor const *, SurfaceTriangle const *>::const_iterator;
//	size_t n = bestPairs.size();
//
//#pragma omp parallel for
//	for (size_t i = 0; i < n; ++i) {
//
//		CompactPatchDescriptor const * cpd1_pt = bestPairs[i].cpd1; /**< ID of the first CPD of the current CPDpair */
//		CompactPatchDescriptor const * cpd2_pt = bestPairs[i].cpd2; /**< ID of the second CPD of the current CPDpair */
//
//		/**
//		 * The equal_range method returns the bounds of a range that includes all the
//		 * elements in the container with a key that compares equal to the given one.
//		 * If the given key does not match any key in the container, the range returned
//		 * has end as both its lower and upper range bounds. The lower and upper bounds
//		 * of the range are returned as a pair of iterators.
//		 */
//		pair<m_it,m_it> range1, range2;
//		range1 = surfTriangleMap_r.equal_range(cpd1_pt);
//		range2 = surfTriangleMap_l.equal_range(cpd2_pt);
//
//		/**
//		 * There is a chance that there might not be any compatible SurfaceTriangles where
//		 * one of the two (or both) CompactPatchDescriptors appear. In this case, just go on
//		 * to the next best matching CPDpair.
//		 */
//		if (range1.first == surfTriangleMap_r.end() && range1.second == surfTriangleMap_r.end())
//			continue;
//		if (range2.first == surfTriangleMap_l.end() && range2.second == surfTriangleMap_l.end())
//			continue;
//
//		for (m_it it1 = range1.first; it1 != range1.second; ++it1) {
//			for (m_it it2 = range2.first; it2 != range2.second; ++it2) {
//				/**
//				 * Sanity-check: the key of all returned values must be equal to the search
//				 * key (cpd1 and cpd2).
//				 */
//				assert(cpd1_pt == it1->first && cpd2_pt == it2->first);
//
//				SurfaceTriangle const * tr1_pt = it1->second; /**< current SurfaceTriangle where the CPD cpd1 appears */
//				SurfaceTriangle const * tr2_pt = it2->second; /**< current SurfaceTriangle where the CPD cpd2 appears */
//				/**
//				 * A SurfaceTriangle has three vertices determined by the centers of
//				 * the three SurfacePatches it points to. Let's call A the center of
//				 * p1, B the center of p2, and C the center of p3. Now, we already
//				 * know that the patch with ID firstPatch_id in the first
//				 * SurfaceTriangle has a good match in terms of Zernike Descriptors
//				 * with the patch with id secondPatch_id in the second SurfaceTriangle.
//				 * Thus, the number of possible surface triangles drops to two, as we
//				 * have already matched two vertices.
//				 */
//				SurfaceTriangle tr1  = tr1_pt->trPermutation1(cpd1_pt);
//				SurfaceTriangle tr21 = tr2_pt->trPermutation1(cpd2_pt);
//				SurfaceTriangle tr22 = tr2_pt->trPermutation2(cpd2_pt);
//
//				DockingPose pose1(tr1, tr21);
//				double score1 = tr1.getSimilarityScore(tr21);
//
//				DockingPose pose2(tr1, tr22);
//				double score2 = tr1.getSimilarityScore(tr22);
//
//
////				point3D nA1  = tr1.p1->patchNormal;
////				point3D nA21 = tr21.p1->patchNormal;
////				point3D nA22 = tr22.p1->patchNormal;
//
//				point3D nA21_T = nA21.transform(pose1.rotation);
//				point3D nA22_T = nA22.transform(pose2.rotation);
//
//				if (nA1.angle(nA21_T) > M_PI/2) {
//#pragma omp critical
//					addPose(pose1, score1, max_num_poses, dockingPoses_map, orderedDockingPoses);
//				}
//				else if(nA1.angle(nA22_T) > M_PI/2) {
//#pragma omp critical
//					addPose(pose2, score2, max_num_poses, dockingPoses_map, orderedDockingPoses);
//				}
//			}
//		}
//	}
//	/*
//	 * Create the output data-structures
//	 */
//	for (auto const & it : dockingPoses_map) {
//		dockingPoses.push_back(it.first);
//		dockingScores.insert(DockingScore(dockingPoses.back(), it.second));
//	}
//};
