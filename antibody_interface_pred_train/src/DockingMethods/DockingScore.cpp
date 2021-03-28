/*
 * DockingScore.cpp
 *
 *  Created on: 25/feb/2015
 *      Author: sebastian
 */


#include "DockingScore.h"

DockingScore::DockingScore() :
		pose(NULL), score(-1) {
}

DockingScore::DockingScore(DockingScore const & ds) :
		pose(ds.pose), score(ds.score) {
}

DockingScore::DockingScore(DockingPose const & p, double score) :
		pose(&p), score(score) {
}

bool DockingScore::operator==(DockingScore const & rhs) const {
	return (doubleCompare(this->score, rhs.score) && *(this->pose) == *(rhs.pose));
}

bool DockingScore::operator<(DockingScore const & rhs) const {
	return (this->score < rhs.score);
}
