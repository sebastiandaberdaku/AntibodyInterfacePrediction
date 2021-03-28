/*
 * stoc.h
 *
 *  Created on: 14/dic/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_STOC_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_STOC_H_

#include <string>

using namespace std;

/**
 * string to char converter.
 * Converts
 */

static inline char stoc(string const & str) {
	if (str.empty())
		return '\0';
	else
		return str[0];
}




#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_STOC_H_ */
