/*
 * node.h
 *
 *  Created on: Sep 24, 2015
 *      Author: sebastian
 */

#ifndef NODE_H_
#define NODE_H_

// The data structures of the proposed improved algorithm.
// The definition of data structure, for improved 2D algorithm.
typedef struct node {
	node();
	node(node const & n);
	node(int zl, int zr, int y, int direction_y);
	node & operator=(node const & n);
	int zl, zr; // zl/zr – the ID tag of the leftmost/rightmost seed of current scan range
	int y; // y–the order number of current scan line
	int direction_y; //the search direction of y axis, -1 -> backward, + 1 -> forward
} node;



#endif /* NODE_H_ */
