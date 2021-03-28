/*
 * sliceNode.h
 *
 *  Created on: Sep 24, 2015
 *      Author: sebastian
 */

#ifndef SLICENODE_H_
#define SLICENODE_H_

#include <vector>

#include "../SeedFill3D/range.h"

using namespace std;


typedef vector<range> IdList;

//The definition of data structure, for proposed 3D algorithm.
typedef struct sliceNode {
	sliceNode();
	sliceNode(sliceNode const & sn);
	sliceNode(IdList const & seedList, int x, int direction_x);
	/**
	 * Copy assignment operator
	 */
	sliceNode & operator=(sliceNode const & sn);
	int x; // each slice is identified by the x coordinate
	IdList seedList; //the seed list of current seed-filling region
	int direction_x; //the search direction of x axis, -1 -> downward, +1 -> upward
} sliceNode;


#endif /* SLICENODE_H_ */
