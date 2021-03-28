/*
 * partialEDMap.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_PARTIALEDMAP_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_PARTIALEDMAP_H_

#include "voxel.h"
#include <map>

/**
 * Struct defining a distanceMapValue data type.
 */
#if defined TIGHT_PACKING
#pragma pack(push)  /* push current alignment to stack */
#pragma pack(1)     /* set alignment to 1 byte boundary */
#endif
typedef struct distanceMapValue {
	distanceMapValue();
	distanceMapValue(distanceMapValue const & value);
	distanceMapValue(uint16_t sDistance, voxel nearestSurfVox);
	distanceMapValue & operator=(distanceMapValue const & value);
	uint16_t sDistance; /**< Squared distance from the nearest boundary/surface voxel. */
	voxel nearestSurfVox; /**< The nearest surface voxel. */
} distanceMapValue;
#if defined TIGHT_PACKING
#pragma pack(pop)   /* restore original alignment from stack */
#endif
/**
 * Struct defining a partial Euclidean Distance Map data type.
 * This data type serves as a container for a map data type.
 */
typedef struct partialEDMap {
	void insert(voxel const & key, distanceMapValue const & value);
	distanceMapValue* find(voxel const & key);
	struct cmp_voxel {
		bool operator()(voxel const & a, voxel const & b);
	};
	std::map<voxel, distanceMapValue, cmp_voxel> distanceMap;
} partialEDMap;

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_PARTIALEDMAP_H_ */
