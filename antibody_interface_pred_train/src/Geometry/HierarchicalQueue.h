/*
 * HierarchicalQueue.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_HIERARCHICALQUEUE_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_HIERARCHICALQUEUE_H_

#include "voxel.h"
#include <queue>

/**
 * Hierarchical queues (HQ) are extremely efficient structures for fast computation
 * of Distance Transforms.
 * A hierarchical queue is the assembly of N simple queues. A queue is also called
 * FIFO ("First In First Out") register. In each queue, tokens are extracted in the
 * same chronological order as they were introduced. Each queue has a single priority
 * level: queue 0 has the highest priority.
 */
typedef struct HierarchicalQueue {
public:
	HierarchicalQueue(uint16_t numberOfQueues);
	virtual ~HierarchicalQueue();
	bool empty();
	size_t size();
	void pop();
	void push(voxel const & vox, uint16_t idx);
	voxel& front();
	const voxel& front() const;
private:
	uint16_t numberOfQueues;
	std::queue<voxel>* queues;
	const std::queue<voxel>& operator[] (size_t idx) const;
	std::queue<voxel>& operator[] (size_t idx);
} HierarchicalQueue;

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_HIERARCHICALQUEUE_H_ */
