/*
 * DockingMethods.h
 *
 *  Created on: Feb 26, 2015
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_DOCKINGMETHODS_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_DOCKINGMETHODS_H_

#include "CPDpair.h"
#include "DockingPose.h"
#include "DockingScore.h"
#include "orderedPosesList.h"
#include "SurfaceTriangle.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>

typedef std::pair<DockingPose, double> scoredDockingPose;

typedef struct scoredDockingPose_compare {
    bool operator() (scoredDockingPose const & lhs, scoredDockingPose const & rhs) const{
        return lhs.second < rhs.second;
    }
} scoredDockingPose_compare;

namespace std {
template<>
struct hash<scoredDockingPose> {
	size_t operator()(scoredDockingPose const & s) const {
		size_t result = hash<DockingPose>()(s.first);
		boost::hash_combine(result, s.second);
		return result;
	}
};
}

typedef struct DockingMethods {
public:
	static void calculateBestComplementarityPoses(vector<CPDpair> const & bestPairs,
			unordered_multimap<CompactPatchDescriptor const *, SurfaceTriangle const *> & surfTriangleMap_r,
			unordered_multimap<CompactPatchDescriptor const *, SurfaceTriangle const *> & surfTriangleMap_l,
			size_t max_num_poses, vector<DockingPose> & dockingPoses_map, orderedPosesList & dockingScores);

private:
	static void addPose(DockingPose const & currentPose, double currentScore, size_t max_num_poses,
			unordered_map<DockingPose, double> & dockingPoses_map,
			set<scoredDockingPose, scoredDockingPose_compare> & orderedDockingPoses);

} DockingMethods;

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_DOCKINGMETHODS_H_ */
