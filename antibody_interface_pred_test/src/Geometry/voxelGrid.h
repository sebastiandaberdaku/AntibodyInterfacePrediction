/*
 * voxelGrid.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_VOXELGRID_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_VOXELGRID_H_

#include "voxel.h"

/**
 * Determine the word size from the preprocessor defines.
 */
#if defined WORD32
/** 32-bit CPU word */
typedef uint32_t word_t;
#else
/** 64-bit CPU word */
typedef uint64_t word_t;
#endif

/**
 * Struct defining the voxelized representation of the molecule.
 */
typedef struct voxelGrid {
private:
	/**
	 * Returns the block index of the voxel with the given indexes.
	 * \param i		ix coordinate of the desired voxel
	 * \param j		iy coordinate of the desired voxel
	 * \param k		iz coordinate of the desired voxel
	 * \return 		the block index of the desired voxel
	 */
	inline size_t getBlockIndex(uint16_t const & i, uint16_t const & j, uint16_t const & k) {
		size_t idx = (i * width + j) * height + k;
		return idx / wordSize;
	};
	/**
	 * Returns the offset value of the voxel within the block.
	 * \param i		ix coordinate of the desired voxel
	 * \param j		iy coordinate of the desired voxel
	 * \param k		iz coordinate of the desired voxel
	 * \return 		the offset value
	 */
	inline size_t getOffset(uint16_t const & i, uint16_t const & j,
			uint16_t const & k) {
		size_t idx = (i * width + j) * height + k;
		return idx % wordSize;
	};
	static const word_t one;
	static const word_t zeros;
	static const uint16_t wordSize;
public:
	/**
	 * Constructor of the voxelGrid object.
	 *
	 * \param length			The length of the voxel grid.
	 * \param width				The width of the voxel grid.
	 * \param height			The height of the voxel grid.
	 */
	voxelGrid(uint16_t length, uint16_t width, uint16_t height);
	/**
	 * Constructor for a cubic grid.
	 * \param dim	length = width = height = dim.
	 */
	voxelGrid(uint16_t dim);
	/**
	 * Copy constructor of the voxelGrid object. Creates a new copy of
	 * an existing voxelGrid.
	 *
	 * \param grid	The voxelGrid to copy.
	 */
	voxelGrid(voxelGrid const & grid);
	/**
	 * Copy assignment operator
	 */
	voxelGrid & operator=(voxelGrid const & grid);

	voxelGrid & operator|=(voxelGrid const & grid);

	/**
	 * Inverts all voxels in the current voxelGrid
	 * @return a read/write reference to the current voxelGrid
	 */
	inline voxelGrid & invert() {
		for (int i = 0; i < this->size; ++i)
			this->flags[i] = ~(this->flags[i]);
		return *this;
	};
	/**
	 * Bitwise NOT operator
	 * @return a read/write reference to the inverted voxelGrid
	 */
	const voxelGrid operator~() const;

	/**
	 * Method for setting the elements of the grid by their 3D coordinates.
	 * \param i		The x coordinate of the voxel to be accessed.
	 * \param j		The y coordinate of the voxel to be accessed.
	 * \param k		The z coordinate of the voxel to be accessed.
	 */
	inline void setVoxel(int32_t const & i, int32_t const & j, int32_t const & k) {
	#if defined RANGECHECK_TEST
		if (i < 0 || i >= length)
			throw std::out_of_range(
					"voxelGrid::setVoxel() - Index i is out of bounds!");
		if (j < 0 || j >= width)
			throw std::out_of_range(
					"voxelGrid::setVoxel() - Index j is out of bounds!");
		if (k < 0 || k >= height)
			throw std::out_of_range(
					"voxelGrid::setVoxel() - Index k is out of bounds!");
	#endif
		size_t idx = (i * width + j) * height + k;
		flags[idx / wordSize] |= (one << (idx % wordSize));
	};
	inline void setVoxel(voxel const & v) {
		size_t idx = (v.ix * width + v.iy) * height + v.iz;
		flags[idx / wordSize] |= (one << (idx % wordSize));
	};
	/**
	 * Method for clearing the elements of the grid by their 3D coordinates.
	 * \param i		The x coordinate of the voxel to be accessed.
	 * \param j		The y coordinate of the voxel to be accessed.
	 * \param k		The z coordinate of the voxel to be accessed.
	 */
	inline void clearVoxel(int32_t const & i, int32_t const & j,
			int32_t const & k) {
	#if defined RANGECHECK_TEST
		if (i < 0 || i >= length)
			throw std::out_of_range(
					"voxelGrid::clearVoxel() - Index i is out of bounds!");
		if (j < 0 || j >= width)
			throw std::out_of_range(
					"voxelGrid::clearVoxel() - Index j is out of bounds!");
		if (k < 0 || k >= height)
			throw std::out_of_range(
					"voxelGrid::clearVoxel() - Index k is out of bounds!");
	#endif
		size_t idx = (i * width + j) * height + k;
		flags[idx / wordSize] &= ~(one << (idx % wordSize));
	};
	inline void clearVoxel(voxel const & v) {
		size_t idx = (v.ix * width + v.iy) * height + v.iz;
		flags[idx / wordSize] &= ~(one << (idx % wordSize));
	};
	/**
	 * Method for reading the elements of the grid by their 3D coordinates.
	 * \param i		The x coordinate of the voxel to be accessed.
	 * \param j		The y coordinate of the voxel to be accessed.
	 * \param k		The z coordinate of the voxel to be accessed.
	 */
	inline bool getVoxel(int32_t const & i, int32_t const & j, int32_t const & k) const {
	#if defined RANGECHECK_TEST
		if (i < 0 || i >= length)
			throw std::out_of_range(
					"voxelGrid::getVoxel() - Index i is out of bounds!");
		if (j < 0 || j >= width)
			throw std::out_of_range(
					"voxelGrid::getVoxel() - Index j is out of bounds!");
		if (k < 0 || k >= height)
			throw std::out_of_range(
					"voxelGrid::getVoxel() - Index k is out of bounds!");
	#endif
		size_t idx = (i * width + j) * height + k;
		return (flags[idx / wordSize] & (one << (idx % wordSize))) != zeros;
	};
	inline bool getVoxel(voxel const & v) const {
		size_t idx = (v.ix * width + v.iy) * height + v.iz;
		return (flags[idx / wordSize] & (one << (idx % wordSize)))
				!= zeros;
	};
	/** Sets all the elements of the grid to false. */
	void clear_all();
	/** Sets all the elements of the grid to true. */
	void set_all();

	bool floodFill(voxel const & startingVoxel);

	friend std::ostream& operator <<(std::ostream& s, const voxelGrid& grid);
	/**
	 * Destructor of the voxelGrid object.
	 */
	virtual ~voxelGrid();
	uint16_t length, width, height;
	size_t size;
	word_t *flags;	/**< Each voxel is composed of a boolean flag. */
} voxelGrid;



#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_VOXELGRID_H_ */
