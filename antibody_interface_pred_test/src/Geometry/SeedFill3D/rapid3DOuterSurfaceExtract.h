/*
 * rapid3DOuterSurfaceExtract.h
 *
 *  Created on: Sep 26, 2015
 *      Author: sebastian
 */

#ifndef RAPID3DOUTERSURFACEEXTRACT_H_
#define RAPID3DOUTERSURFACEEXTRACT_H_

#include "../voxelGrid.h"
#include <stack>
#include "../SeedFill3D/node.h"
#include "../SeedFill3D/sliceNode.h"

namespace rapid3DOuterSurfaceExtract {
/**
 * Clears a whole range in temp, from zl to zr, on line y of slice z.
 * @param cpkModel	the volumetric model
 * @param temp		a temporary data structure
 * @param surface	output voxelGrid
 * @param zl		range lower bound
 * @param zr		range upper bound
 * @param y			range line
 * @param z			image slice
 */
static inline void clearRange(voxelGrid const & cpkModel, voxelGrid & temp, voxelGrid & surface, int x, int y, int zl, int zr) {
	for (int i = zl; i <= zr; i++) {
		if (x < cpkModel.length - 1 && cpkModel.getVoxel(x + 1, y, i))
			surface.setVoxel(x + 1, y, i);
		if (x > 0 && cpkModel.getVoxel(x - 1, y, i))
			surface.setVoxel(x - 1, y, i);
		if (y < cpkModel.width - 1 && cpkModel.getVoxel(x, y + 1, i))
			surface.setVoxel(x, y + 1, i);
		if (y > 0 && cpkModel.getVoxel(x, y - 1, i))
			surface.setVoxel(x, y - 1, i);
		temp.clearVoxel(x, y, i);
	}
	if (zr < surface.height - 1)
		surface.setVoxel(x, y, zr + 1);
	if (zl > 0)
		surface.setVoxel(x, y, zl - 1);
};

/**
 * Returns the unfilled and valid range discovered scanning the target image from seed
 * @param cpkModel		input 3D image
 * @param surface		surface 3D image
 * @param seed		the seed
 * OUTPUT:
 * @param zl 		the left range delimiter
 * @param zr 		the right range delimiter
 */
static inline void obtainUnfilledValidRange(voxelGrid const & cpkModel, voxelGrid const & temp, voxel const & seed,
		int & zl, int & zr) {
	zr = seed.iz;
	while ((zr < cpkModel.height - 1) && !cpkModel.getVoxel(seed.ix, seed.iy, zr + 1) && temp.getVoxel(seed.ix, seed.iy, zr + 1)) {
		++zr;
	}
	zl = seed.iz;
	while ((zl > 0) && !cpkModel.getVoxel(seed.ix, seed.iy, zl - 1) && temp.getVoxel(seed.ix, seed.iy, zl - 1)) {
		--zl;
	}
};

/**
 * Checks if the current range [zpl, zpr] contains any valid seeds. If a seed is found
 * the corresponding range is extracted, filled and inserted in dList.
 * @param cpkModel	the current voxelGrid
 * @param surface	output voxelGrid
 * @param zpl		lower bound of the range to be checked
 * @param zpr		upper bound of the range to be checked
 * @param y			line of the current range
 * @param z			slice of the current range
 * OUTPUT:
 * @param dList		list of the filled ranges in the current slice
 * @param zl		lower bound of the newly found range
 * @param zr		upper bound of the newly found range
 * @return			true if a valid seed is found, false otherwise
 */
static inline bool checkAndFillRange(voxelGrid const & cpkModel, voxelGrid & temp, voxelGrid & surface, int zpl, int zpr, int y, int x,
		IdList & dList, int & zl, int & zr) {
	for (int i = zpl; i <= zpr; ++i) {
		if (!cpkModel.getVoxel(x, y, i) && temp.getVoxel(x, y, i)) { // if a valid seed is found
			// extract the unfilled range,and write the filled seeds ID into dList
			// search the valid and unfilled range ([zl,zr]) from [zpl, zpr]
			obtainUnfilledValidRange(cpkModel, temp, voxel(x, y, i), zl, zr);
			//extract and fill the range [zl,zr], mark the seeds and write the ID of them into dList
			clearRange(cpkModel, temp, surface, x, y, zl, zr);
			dList.push_back(range(zl, zr, y));
			return true;
		}
	}
	return false;
};

/**
 * This method flood-fills a single slice of the input voxelGrid, starting from voxel seed.
 * The slice is identified by the z coordinate of seed. The ID of each seed is extracted from
 * dList.
 * @param cpkModel	target voxelGrid
 * @param surface	output voxelGrid
 * @param seed		starting voxel for the flood-filling procedure
 * OUTPUT:
 * @param seedList	list of the filled ranges in the current slice
 */
static inline void seedFill2D(voxelGrid const & cpkModel, voxelGrid & temp, voxelGrid & surface, voxel const & seed, IdList & seedList) {
	// empty the seed list
	seedList.clear();
	// initialize the empty 2D stack
	stack<node> stack2D;
	// obtain the unfilled and valid range ([zl,zr]) by using seed
	int zl, zr; // range [zl, zr]
	obtainUnfilledValidRange(cpkModel, temp, seed, zl, zr);
	// extract and fill the range [zl,zr],
	clearRange(cpkModel, temp, surface, seed.ix, seed.iy, zl, zr);
	// mark the seeds and write the ID of them into dList
	seedList.push_back(range(zl, zr, seed.iy));
	// create two nodes with range ([zl, zr]), push them into 2D stack, along two opposite directions
	int sy = seed.iy, max_width = cpkModel.width;
	if (sy > 0)
		stack2D.push(node(zl, zr, seed.iy, -1));
	if (sy < max_width - 1)
		stack2D.push(node(zl, zr, seed.iy, 1));
	// obtain the length of this range
	while (!stack2D.empty()) { // if the 2D stack is not empty
		//pop a node and then set the search range ([zpl,zpr]) onto the next scan line
		node c(stack2D.top()); // current node
		stack2D.pop();
		int zpl = c.zl, zpr = c.zr, y = c.y + c.direction_y;

		while (checkAndFillRange(cpkModel, temp, surface, zpl, zpr, y, seed.ix, seedList, zl, zr)) { //if the range [zpl, zpr] is valid
			if (zr < zpr - 1) { //it maybe exist other valid and unfilled ranges on the same scan line
				for (int i = zr + 2; i <= zpr; ++i) { // search the leftmost seed of the next valid range, which locates at the same scan line
					if (!cpkModel.getVoxel(seed.ix, y, i) && temp.getVoxel(seed.ix, y, i)) { /* a valid seed exists */
						stack2D.push(node(zr + 2, zpr, y - c.direction_y, c.direction_y));
						break;
					}
				}
			}
			// execute the necessary rollback operation
			if (zl < zpl - 1) {	//the rollback operation may occur on the left side of zpl
			// rollback to the previous scan line
				int old_y = y - c.direction_y;
				// search a valid and unfilled seed on the left side of zpl
				for (int i = zl; i <= zpl - 2; ++i) {
					if (!cpkModel.getVoxel(seed.ix, old_y, i) && temp.getVoxel(seed.ix, old_y, i)) { /* a valid seed exists */
						//the search range is valid, push it into 2D stack
						// push a new node into 2D stack, the range is [zl, zpl-2], the direction is opposite
						stack2D.push(node(zl, zpl - 2, y, -c.direction_y));
						break;
					}
				}
			}
			if (zr > zpr + 1) {	//the rollback operation may occur on the right side of zpr
			// rollback to the previous scan line
				int old_y = y - c.direction_y;
				// search a valid and unfilled seed on the right side of zpr
				for (int i = zpr + 2; i <= zr; ++i) {
					// if (a valid seed exists){//the search range is valid, push it into 2D stack
					if (!cpkModel.getVoxel(seed.ix, old_y, i) && temp.getVoxel(seed.ix, old_y, i)) { /* a valid seed exists */
						// the search range is valid, push it into 2D stack
						// push a new node into 2D stack, the range is [zpr+ 2, zr], the direction is opposite
						stack2D.push(node(zpr + 2, zr, y, -c.direction_y));
						break;
					}
				}
			}
			// continue the loop from this range [zl, zr], along the same direction as before
			// obtain the new search range by assigning the range [zl, zr] to [zpl, zpr]
			zpl = zl;
			zpr = zr;
			y += c.direction_y;
			if (y < 0 || y >= cpkModel.width)
				break;
		}
	}
};

/**
 * This method projects the seed list of the current slice onto slice z and returns true if
 * a valid seed is found in the projected area, false otherwise.
 * @param cpkModel	the current voxelGrid
 * @param temp		temporary voxelGrid
 * @param seed_list	the list of seeds in the current slice
 * @param z			the neighbor slice
 * @param leap_var	the leaping variable is used to leap over plenty of invalid
 * 					voxels in order to accelerate search process for the first
 * 					valid seed
 * OUTPUT:
 * @param sx 		x coordinate of the first valid seed found
 * @param sy		y coordinate of the first valid seed found
 * @return			true if a new seed is found, false otherwise
 */
static inline bool searchNeighborSlice(voxelGrid const & cpkModel, voxelGrid const & temp, IdList & seed_list,
		int x, int leap_var, int & sy, int & sz) {
	IdList::iterator id, begin = seed_list.begin(),  end = seed_list.end();
	for (id = begin; id != end; ++id) {
		for (int iz = id->zl; iz <= id->zr; iz += leap_var) {
			if (!cpkModel.getVoxel(x, id->y, iz) && temp.getVoxel(x, id->y, iz)) {
				sz = iz;
				sy = id->y;
				/*
				 * Resize the current seed, as all voxels from zl to the current ix are occupied
				 */
				id->zl = iz;
				/*
				 * Erase previous seeds, excluding the current one!
				 */
				seed_list.erase(begin, id);
				return true;
			}
		}
	}
	seed_list.clear();
	return false;
};

/**
 * Here we implement a rapid 3D seed-filling algorithm, to extract or fill the object-connected 3D region.
 * The algorithm uses an improved 2D seed-filling algorithm, which extracts connected region in slice quickly and
 * consumes fewer stack operations and less memory compared with the existing algorithms. The improved 2D algorithm
 * is enclosed as a basic unit within the framework of the proposed 3D seed-filling algorithm, in order to reduce
 * the complexity of direction of seeds search, and accelerate region search on adjacent slices. Experimental results
 * demonstrate advantages of this algorithm including eliminating the redundancy of seeds search, repetition of stack
 * operations and running with high efficiency.
 * @param cpkModel input image
 * @param surface	output voxelGrid
 * @param seed		starting voxel
 * @param leap_var	the leaping variable is used to leap over plenty of invalid
 * 					voxels in order to accelerate search process for the first
 * 					valid seed
 */
static inline void outerSurfaceExtract(voxelGrid const & cpkModel, voxelGrid & surface, voxel const & seed, int leap_var = 1) {
	if (cpkModel.getVoxel(seed))
		return;
	voxelGrid temp(cpkModel.length, cpkModel.width, cpkModel.height);
	temp.set_all();
	//create empty 3D stack
	stack<sliceNode> stack3D;

	IdList currentSeeds;
	//extract and fill the initial unfilled region by using an initial seed (seed)
	seedFill2D(cpkModel, temp,  surface, seed, currentSeeds);
	//create two nodes with this initial region(dList), push them into 3D stack, along two opposite directions
	if (seed.ix < cpkModel.length - 1)
		stack3D.push(sliceNode(currentSeeds, seed.ix, +1));
	if (seed.ix > 0)
		stack3D.push(sliceNode(currentSeeds, seed.ix, -1));
	while (!stack3D.empty()) {
		//pop a node from 3D stack
		sliceNode n(stack3D.top());
		stack3D.pop();
		//set the search region (PassList) by assigning the region (dList) to it
		currentSeeds = n.seedList; /*
							 * std::list::operator= c++11 Assigns new contents to the container,
							 * replacing its current contents, and modifying its size accordingly.
							 */

		int new_x = n.x + n.direction_x;
		//obtain a search region (Plist_mapped) by mapping the region (PassList) onto the next slice
		//and then check the validation of this mapped region (Plist_mapped)
		int sy, sz; // new seed coordinates
		while(searchNeighborSlice(cpkModel, temp, currentSeeds, new_x, leap_var, sy, sz)){ //if the region is valid
			IdList newSeeds;
			//search a valid seed in the mapped region (Plist_mapped)
			//extract and fill this new unfilled region (newSeeds here represents the first new region)
			seedFill2D(cpkModel, temp, surface, voxel(new_x, sy, sz), newSeeds);
			//continue to search other new unfilled regions in the region (Plist_mapped)
			while (searchNeighborSlice(cpkModel, temp, currentSeeds, new_x, leap_var, sy, sz)){
				//extract and fill these new unfilled regions
				IdList otherSeeds;
				seedFill2D(cpkModel, temp, surface, voxel(new_x, sy, sz), otherSeeds);
				//create two nodes with each region, push them into 3D stack, along two opposite directions
				if (new_x < cpkModel.length - 1)
					stack3D.push(sliceNode(otherSeeds, new_x, +1));
				if (new_x > 0)
					stack3D.push(sliceNode(otherSeeds, new_x, -1));
			}
			//execute the rollback operation
			//map the region (CurrentList) onto the previous slice
			currentSeeds = newSeeds;
			while (searchNeighborSlice(cpkModel, temp, newSeeds, new_x - n.direction_x, leap_var, sy, sz)) {
				IdList otherSeeds;
				//extract and fill these new unfilled regions
				seedFill2D(cpkModel, temp, surface, voxel(new_x - n.direction_x, sy, sz), otherSeeds);
				//create two nodes with each region, push them into 3D stack, along two opposite directions
				if (new_x - n.direction_x < cpkModel.length - 1)
					stack3D.push(sliceNode(otherSeeds, new_x - n.direction_x, +1));
				if (new_x - n.direction_x > 0)
					stack3D.push(sliceNode(otherSeeds, new_x - n.direction_x, -1));
			}
			//continue the loop from the region (CurrentList), along the same direction as before
			//obtain the new search region by assigning the region (CurrentList) to (PassList)
			new_x += n.direction_x;
			if (new_x < 0 || new_x == cpkModel.length)
				break;
		}
	}
};

}



#endif /* RAPID3DOUTERSURFACEEXTRACT_H_ */
