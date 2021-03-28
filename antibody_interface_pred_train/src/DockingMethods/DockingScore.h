/*
 * DockingScore.h
 *
 *  Created on: 13/feb/2015
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_DOCKINGSCORE_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_DOCKINGSCORE_H_
#include "DockingPose.h"
#include <boost/functional/hash.hpp>

typedef struct DockingScore {
	DockingScore();
	DockingScore(DockingScore const & ds);
	DockingScore(DockingPose const & p, double score);
	DockingScore(DockingPose const * p, double score);

	DockingPose const * pose;
	double score;

	bool operator==(DockingScore const & rhs) const;
	bool operator<(DockingScore const & rhs) const;
} DockingScore;

namespace std {
template<>
struct hash<DockingScore> {
	size_t operator()(DockingScore const & s) const {
		size_t result;
		boost::hash_combine(result, s.pose);
		boost::hash_combine(result, s.score);
		return result;
	}
};
}

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_DOCKINGSCORE_H_ */
