/*
 * trim.h
 *
 *  Created on: 07/set/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_TRIM_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_TRIM_H_

#include <string>

/**
 * Trims the leading and ending space characters from the input string.
 * \param str	The string to be trimmed.
 * \return		The resulting trimmed string.
 */
static inline std::string trim(std::string& str) {
	str.erase(0, str.find_first_not_of(' ')); /* remove prefix spaces */
	str.erase(str.find_last_not_of(' ') + 1); /* remove suffix spaces */
	return str;
};

static inline std::string trim(std::string const & str) {
	std::string temp(str);
	temp.erase(0, temp.find_first_not_of(' ')); /* remove prefix spaces */
	temp.erase(temp.find_last_not_of(' ') + 1); /* remove suffix spaces */
	return temp;
};

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_TRIM_H_ */
