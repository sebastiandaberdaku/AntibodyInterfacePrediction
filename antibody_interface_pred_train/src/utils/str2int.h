/*
 * str2int.h
 *
 *  Created on: 01/dic/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_STR2INT_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_STR2INT_H_

inline constexpr unsigned int str2int(const char* str, int h = 0) {
	return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
}

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_STR2INT_H_ */
