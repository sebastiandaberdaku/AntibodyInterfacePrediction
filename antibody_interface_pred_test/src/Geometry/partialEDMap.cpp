/*
 * partialEDMap.cpp
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

/**
 * Implementation of distanceMapValue and partialEDMap.
 */
#include "partialEDMap.h"

/**
 * Default constructor of distanceMapValue.
 */
distanceMapValue::distanceMapValue() :
		sDistance(0), nearestSurfVox(voxel(0, 0, 0)) { }
/**
 * Copy constructor of distanceMapValue.
 */
distanceMapValue::distanceMapValue(distanceMapValue const & value) :
		sDistance(value.sDistance), nearestSurfVox(value.nearestSurfVox) { }
/**
 * Constructor of distanceMapValue.
 */
distanceMapValue::distanceMapValue(uint16_t sDistance, voxel nearestSurfVox) :
		sDistance(sDistance), nearestSurfVox(nearestSurfVox) { }
/**
 * Copy assignment operator
 */
distanceMapValue & distanceMapValue::operator=(const distanceMapValue & value) {
	if (this != &value) {
		this->nearestSurfVox = value.nearestSurfVox;
		this->sDistance = value.sDistance;
	}
	return *this;
}
/**
 * This method inserts the input <key, value> couple into the
 * distance map. If another couple with the same key exists,
 * the method replaces it.
 * \param key		The voxel corresponding to the current distanceMapValue
 * \param value		The distanceMapValue for the current voxel.
 */
void partialEDMap::insert(voxel const & key, distanceMapValue const & value) {
	std::map<voxel, distanceMapValue, cmp_voxel>::iterator it;
	/* Returns an iterator to the couple with the given key if it exists,
	 * or to map.end() if no existing correspondence is found. */
	it = distanceMap.find(key);
	/* Remove the existing entry if its distance
	 * value is greater than the current one. */
	if (it != distanceMap.end()) {
		distanceMap.erase(it);
	}
	distanceMap.insert(std::pair<voxel, distanceMapValue>(key, value));
}
/**
 * Simple find method. Returns a pointer to the corresponding distance
 * map value if an entry with the given key already exists in the map,
 * NULL otherwise.
 * \param key 	The key of the desired distanceMapValue.
 * \return		A reference to the desired distanceMapValue if it
 * 				exists, NULL otherwise.
 */
distanceMapValue* partialEDMap::find(voxel const & key) {
	std::map<voxel, distanceMapValue, cmp_voxel>::iterator it;
	it = distanceMap.find(key);
	if (it == distanceMap.end())
		return NULL;
	return &(it->second);
}
/**
 * Comparator for the voxel data type.
 */
bool partialEDMap::cmp_voxel::operator()(voxel const & a, voxel const & b) {
	if (a.ix != b.ix)
		return (a.ix < b.ix);
	else if (a.iy != b.iy)
		return (a.iy < b.iy);
	else
		return (a.iz < b.iz);
}
