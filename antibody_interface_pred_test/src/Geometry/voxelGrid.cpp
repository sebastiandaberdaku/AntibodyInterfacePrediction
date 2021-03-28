/*
 * voxelGrid.cpp
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#include "voxelGrid.h"
#include <parallel/algorithm>
#include <queue>

const uint16_t voxelGrid::wordSize = 8 * sizeof(word_t);
const word_t voxelGrid::one = 1;
const word_t voxelGrid::zeros = 0;

/**
 * Constructor of the voxelGrid object.
 *
 * \param length			The length of the voxel grid.
 * \param width				The width of the voxel grid.
 * \param height			The height of the voxel grid.
 */
voxelGrid::voxelGrid(uint16_t length, uint16_t width, uint16_t height) :
		length(length), width(width), height(height) {
	/* size in bit */
	double size_b = length * width * height;
	size = static_cast<size_t>(ceil(size_b / wordSize));
	flags = new word_t[size];
	clear_all();
}
/**
 * Constructor for a cubic grid.
 * \param dim	length = width = height = dim.
 */
voxelGrid::voxelGrid(uint16_t dim) :
		length(dim), width(dim), height(dim) {
	/* size in bit */
	double size_b = length * width * height;
	size = static_cast<size_t>(ceil(size_b / wordSize));
	flags = new word_t[size];
	clear_all();
}
/**
 * Copy constructor of the voxelGrid object. Creates a new copy of
 * an existing voxelGrid.
 *
 * \param grid	The voxelGrid to copy.
 */
voxelGrid::voxelGrid(voxelGrid const & grid) :
		length(grid.length), width(grid.width), height(grid.height), size(grid.size) {
	flags = new word_t[size]();
	std::copy(grid.flags, grid.flags + grid.size, this->flags);
}
/**
 * Copy assignment operator
 */
voxelGrid & voxelGrid::operator=(voxelGrid const & grid) {
	if (this != &grid) {// protect against invalid self-assignment
		this->length = grid.length;
		this->width = grid.width;
		this->height = grid.height;
		this->size = grid.size;
        // 1: allocate new memory and copy the elements
		word_t * new_flags = new word_t[this->size]();
		std::copy(grid.flags, grid.flags + grid.size, new_flags);
        // 2: deallocate old memory
		if (this->flags != NULL)
			delete[] this->flags;
        // 3: assign the new memory to the object
		this->flags = new_flags;
	}
    // by convention, always return *this
	return *this;
}
voxelGrid & voxelGrid::operator|=(voxelGrid const & grid) {
	if (grid.size != this->size)
		return *this;
	for (int i = 0; i < this->size; ++i)
		this->flags[i] = this->flags[i] | grid.flags[i];
	return *this;
}

/**
 * Bitwise NOT operator
 * @return a read/write reference to the inverted voxelGrid
 */
const voxelGrid voxelGrid::operator~() const {
	return (voxelGrid(*this).invert());
}
/** Sets all the elements of the grid to false. */
void voxelGrid::clear_all() {
	std::fill(flags, flags + size, zeros);
}

/** Sets all the elements of the grid to true. */
void voxelGrid::set_all() {
	std::fill(flags, flags + size, ~zeros);
}

bool voxelGrid::floodFill(voxel const & startingVoxel) {
	if (this->getVoxel(startingVoxel))
		return false;
	queue<voxel> pendingVoxels;
	pendingVoxels.push(startingVoxel);

	voxel c; // current voxel

	while (!pendingVoxels.empty()) {
		c = pendingVoxels.front();
		pendingVoxels.pop();
		for (int ii = c.ix - 1; ii >= 0 && !this->getVoxel(ii, c.iy, c.iz); --ii) {
			this->setVoxel(ii, c.iy, c.iz);
			if ((c.iy - 1 > 0) && !this->getVoxel(ii, c.iy - 1, c.iz))
				pendingVoxels.push(voxel(ii, c.iy - 1, c.iz));
			if ((c.iy + 1 < this->width) && !this->getVoxel(ii, c.iy + 1, c.iz))
				pendingVoxels.push(voxel(ii, c.iy + 1, c.iz));
			if ((c.iz - 1 > 0) && !this->getVoxel(ii, c.iy, c.iz - 1))
				pendingVoxels.push(voxel(ii, c.iy, c.iz - 1));
			if ((c.iz + 1 < this->height) && !this->getVoxel(ii, c.iy, c.iz + 1))
				pendingVoxels.push(voxel(ii, c.iy, c.iz + 1));
		}
		for (int ii = c.ix; ii < this->length && !this->getVoxel(ii, c.iy, c.iz); ++ii) {
			this->setVoxel(ii, c.iy, c.iz);
			if ((c.iy - 1 > 0) && !this->getVoxel(ii, c.iy - 1, c.iz))
				pendingVoxels.push(voxel(ii, c.iy - 1, c.iz));
			if ((c.iy + 1 < this->width) && !this->getVoxel(ii, c.iy + 1, c.iz))
				pendingVoxels.push(voxel(ii, c.iy + 1, c.iz));
			if ((c.iz - 1 > 0) && !this->getVoxel(ii, c.iy, c.iz - 1))
				pendingVoxels.push(voxel(ii, c.iy, c.iz - 1));
			if ((c.iz + 1 < this->height) && !this->getVoxel(ii, c.iy, c.iz + 1))
				pendingVoxels.push(voxel(ii, c.iy, c.iz + 1));
		}
	}
	return true;
}

std::ostream& operator <<(std::ostream& s, const voxelGrid& grid) {
	for (uint64_t i = 0; i < grid.size; ++i) {
		s << grid.flags[i];
	}
	return s;
}
/**
 * Destructor of the voxelGrid object.
 */
voxelGrid::~voxelGrid() {
	delete[] flags;
}

