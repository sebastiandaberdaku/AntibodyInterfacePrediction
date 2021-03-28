/*
 * numerical_limits.h
 *
 *  Created on: 15/dic/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_NUMERICAL_LIMITS_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_NUMERICAL_LIMITS_H_

#include <limits>

const float max_float = std::numeric_limits<float>::max();
const float min_float = -max_float;

const double max_double = std::numeric_limits<double>::max();
const double min_double = -max_double;

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_NUMERICAL_LIMITS_H_ */
