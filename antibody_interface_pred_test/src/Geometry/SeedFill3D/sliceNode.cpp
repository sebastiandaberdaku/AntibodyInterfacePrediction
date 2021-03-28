/*
 * sliceNode.cpp
 *
 *  Created on: Sep 24, 2015
 *      Author: sebastian
 */

#include "sliceNode.h"

sliceNode::sliceNode() : x(0), direction_x(0) { };
sliceNode::sliceNode(sliceNode const & sn) :
		x(sn.x), seedList(sn.seedList), direction_x(sn.direction_x) { };
sliceNode::sliceNode(IdList const & seedList, int x, int direction_x) :
		x(x), seedList(seedList), direction_x(direction_x) { };
/**
 * Copy assignment operator
 */
sliceNode & sliceNode::operator=(sliceNode const & sn) {
	if (this != &sn) { // protect against invalid self-assignment
		this->x = sn.x;
		this->seedList = sn.seedList;
		this->direction_x = sn.direction_x;
	}
	// by convention, always return *this
	return *this;
}

