/*
 * orderedScoreList.h
 *
 *  Created on: Feb 15, 2015
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_ORDEREDPOSESLIST_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_ORDEREDPOSESLIST_H_

#include "DockingScore.h"
#include <vector>
#include <math.h>
#include <numeric>

using namespace std;

class orderedPosesList {
public:
	vector<DockingScore> scores_list;
private:
	/**
	 * Simple binary search algorithm used to keep the list of docking scores ordered.
	 * The method returns the index where the current key must be inserted. The list is
	 * ordered in decreasing order of the score entry of DockingPose. The search is done in
	 * O(log n) time, where n is the length of the input DockingPose vector.
	 *
	 * @param key		DockingPose whose position we want to find
	 * @return			index where the current DockingPose must be inserted in order to keep
	 * 					the scores_list vector ordered
	 */
	inline size_t binarySearch(DockingScore const & key) {
		if (scores_list.empty())
			return 0;
		size_t lo = 0;
		size_t hi = scores_list.size();
		// To start, find the index of the middle position.
		size_t position = (lo + hi) / 2;

		while (!doubleCompare(scores_list[position].score, key.score) && (lo < hi)) {
			if (scores_list[position].score < key.score) {
				hi = position;
			} else {
				lo = position + 1;
			}
			position = (lo + hi) / 2;
		}
		return position;
	};

public:
	size_t size() const;
	bool empty() const;
	void clear();
	void reserve(size_t n);
	void insert(DockingScore const & docking_score);
	double bestScore() const;
	double worstScore() const;
	void kBestScores(int k, vector<DockingScore> & bestScores) const;
	DockingScore & operator[](size_t idx);
	DockingScore const & operator[](size_t idx) const;
	void kBestScores(int k, orderedPosesList & bestScores);
};

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_DOCKINGMETHODS_ORDEREDPOSESLIST_H_ */
