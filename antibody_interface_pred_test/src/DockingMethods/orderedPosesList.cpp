/*
 * orderedPosesList.cpp
 *
 *  Created on: 28/feb/2015
 *      Author: sebastian
 */

#include "orderedPosesList.h"
/**
 * Gives the number of elements in the orderedPosesList
 * @return number of poses in the current list
 */
size_t orderedPosesList::size() const {
	return scores_list.size();
}
/**
 * Returns true if the list is empty, false otherwise.
 * @return	true if empty, false otherwise
 */
bool orderedPosesList::empty() const {
	return scores_list.empty();
}
/**
 * Removes all elements from the current list.
 */
void orderedPosesList::clear() {
	scores_list.clear();
}
/**
 * Requests that the capacity be at least enough to contain n elements.
 * If n is greater than the current vector capacity, the function causes
 * the container to reallocate its storage increasing its capacity to n
 * (or greater). In all other cases, the function call does not cause a
 * reallocation and the vector capacity is not affected. This function
 * has no effect on the vector size and cannot alter its elements.
 *
 * @param n 	Minimum capacity for the vector. Note that the resulting
 * 				vector capacity may be equal or greater than n. Member
 * 				type size_type is an unsigned integral type.
 */
void orderedPosesList::reserve(size_t n) {
	scores_list.reserve(n);
}
/**
 * Inserts an element in the list. The method determines the correct
 * position for the element to be inserted in order to keep the list
 * ordered.
 * @param docking_score		element to be inserted
 */
void orderedPosesList::insert(DockingScore const & docking_score) {
	size_t idx = binarySearch(docking_score);
	scores_list.insert(scores_list.begin() + idx, docking_score);
}
/**
 * Returns the best score in the list. If the list is empty the behavior is undefined.
 * @return	best score in the list.
 */
double orderedPosesList::bestScore() const {
	return scores_list.front().score;
}
/**
 * Returns the worst score in the list. If the list is empty the behavior is undefined.
 * @return	worst score in the list.
 */
double orderedPosesList::worstScore() const{
	return scores_list.back().score;
}
/**
 * This method extracts approximately the k-best elements from the current list.
 * Elements in the list are ordered in decreasing order for increasing indexes.
 * The algorithm first extracts elements from index 0 to k-1 and calculates
 * the mean value of the interval and the standard deviation. After that,
 * elements surrounding index k are evaluated. An element with index greater
 * than k-1 can be added to the list if its score is greater (or equal) than
 * mean_value - std_dev. Likewise, an element with index smaller or equal to
 * k-1 can be removed from the list if its score is smaller than
 * mean_value - std_dev. This way, the extraction procedure of the top k elements
 * does not miss good poses near the cut-off index, and does not introduce
 * low-scored poses.
 *
 * @param k 	approximate number of elements to extract
 * @param bestScores	output vector of the k-best scores
 */
void orderedPosesList::kBestScores(int k, vector<DockingScore> & bestScores) const{
	if (scores_list.empty() || k == 0 || scores_list.size() < k)
		return;

	vector<double> v(k);

	for (int i = 0; i < k; ++i)
		v[i] = scores_list[i].score;

	double sum = accumulate(v.begin(), v.end(), (double)0.0);
	double mean = sum / k;

	vector<double> diff(k);
	transform(v.begin(), v.end(), diff.begin(), bind2nd(minus<double>(), mean));
	double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	double std_dev;

	if (k == 1)
		std_dev = 0;
	else
		std_dev = sqrt(sq_sum / (k - 1)); //  corrected sample standard deviation

	double threshold = mean - std_dev; // or mean - 2 * std_dev

	while (k > 0 && scores_list[k - 1].score < threshold)
		--k;
	while (k < scores_list.size() && scores_list[k].score >= threshold)
		++k;

	bestScores.resize(k);
	copy(scores_list.begin(), scores_list.begin() + k, bestScores.begin());
}
/**
 * Subscript operator
 * @param idx	index of the desired element
 * @return	a read/write reference to the element at index idx
 */
DockingScore & orderedPosesList::operator[](size_t idx) {
	return (scores_list[idx]);
}
/**
 * Subscript operator
 * @param idx	index of the desired element
 * @return	a const reference to the element at index idx
 */
DockingScore const & orderedPosesList::operator[](size_t idx) const {
	return (scores_list[idx]);
}

/**
 * This method extracts approximately the k-best elements from the current list.
 * Elements in the list are ordered in decreasing order for increasing indexes.
 * The algorithm first extracts elements from index 0 to k-1 and calculates
 * the mean value of the interval and the standard deviation. After that,
 * elements surrounding index k are evaluated. An element with index greater
 * than k-1 can be added to the list if its score is greater (or equal) than
 * mean_value - std_dev. Likewise, an element with index smaller or equal to
 * k-1 can be removed from the list if its score is smaller than
 * mean_value - std_dev. This way, the extraction procedure of the top k elements
 * does not miss good poses near the cut-off index, and does not introduce
 * low-scored poses.
 *
 * @param k 	approximate number of elements to extract
 * @param bestScores	output orderedPosesList of the k-best scores
 */
void orderedPosesList::kBestScores(int k, orderedPosesList & bestScores) {
	kBestScores(k, bestScores.scores_list);
}


