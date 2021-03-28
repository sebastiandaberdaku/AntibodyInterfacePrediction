/*
 * CPDpair.h
 *
 *  Created on: 20/ago/2014
 *      Author: sebastian
 */

#ifndef CPDPAIR_H_
#define CPDPAIR_H_

#include "../SurfacePatch/CompactPatchDescriptor.h"
#include <vector>
#include <functional>
#include <boost/functional/hash.hpp>

using namespace std;

/**
 * This struct defines a pair of matching CompactPatchDescriptors on the
 * two molecular surfaces we want to dock. It actually contains only pointers
 * to the matching CPDs, and a score variable which gives a measure of the
 * quality of the match.
 */
typedef struct CPDpair {
	/** pointers to the CPDs composing the current pair */
	const CompactPatchDescriptor *cpd1, *cpd2;
	/** score of the current pair's match */
	double score_surface, score_electrostatics, score_combined;

	CPDpair();
	CPDpair(CPDpair const & p);
	CPDpair(CompactPatchDescriptor const & cpd1, CompactPatchDescriptor const & cpd2);

	/**
	 * Copy assignment operator
	 */
	CPDpair & operator=(CPDpair const & p);
	/**
	 * Equal?
	 */
	inline bool operator==(CPDpair const & rhs) const {
		return (*(this->cpd1) == *rhs.cpd1) && (*(this->cpd2) == *rhs.cpd2);
	};

	friend inline bool operator>(CPDpair const & lhs, CPDpair const & rhs) {
		return lhs.score_combined > rhs.score_combined;
	};
	friend inline bool operator>=(CPDpair const & lhs, CPDpair const & rhs) {
		return lhs.score_combined >= rhs.score_combined;
	};
	friend inline bool operator<(CPDpair const & lhs, CPDpair const & rhs) {
		return lhs.score_combined < rhs.score_combined;
	};
	friend inline bool operator<=(CPDpair const & lhs, CPDpair const & rhs) {
		return lhs.score_combined <= rhs.score_combined;
	};

	friend ostream & operator<<(ostream &os, CPDpair const & p);
} CPDpair;

namespace std {
template<>
struct hash<CPDpair> {
	size_t operator()(CPDpair const & p) const {
		size_t result = hash<size_t>()(p.cpd1->ID);
		boost::hash_combine(result, p.cpd2->ID);
		return result;
	}
};

template<>
struct equal_to<CPDpair> {
	bool operator()(CPDpair const & p, CPDpair const & q) const {
		return (*p.cpd1 == *q.cpd1) && (*p.cpd2 == *q.cpd2);
	}
};
}

#endif /* CPDPAIR_H_ */
