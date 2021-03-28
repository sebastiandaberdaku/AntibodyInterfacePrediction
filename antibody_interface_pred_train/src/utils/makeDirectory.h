/*
 * makeDirectory.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_MAKEDIRECTORY_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_MAKEDIRECTORY_H_

#include <boost/regex.hpp>
#include <cerrno>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;
using namespace boost;
/**
 * makeDirectory() attempts to create a directory named pathname.
 * This method uses a simple regex rule to check if the given string is
 * a valid UNIX pathname. Then it checks with stat if the directory already
 * exists. Finally, it tries to create the directory with mkdir.
 *
 * @param pathname	the name of the directory to create
 *
 * @see man 2 stat
 * @see man 2 mkdir
 */
static inline void makeDirectory(char const * pathname) {
	static regex r("[^[.NUL.]]+");
	if (!regex_match(pathname, r))
		throw invalid_argument("makeDirectory::makeDirectory() - The given string is not a valid pathname.");
	// stat for checking if the directory exists,
	// mkdir, to create the directory.
	struct stat st = { 0 };
	if (stat(pathname, &st) == -1) {
		if (mkdir(pathname, 0700) == -1) {
			char const * error_message;
			switch (errno) {
			case EACCES:
				error_message =  "The parent directory does not allow write permission to the process, or one of the directories in pathname did not allow search permission.";
				break;
			case EDQUOT:
				error_message =	"The user's quota of disk blocks or inodes on the filesystem has been exhausted.";
				break;
			case EEXIST:
				error_message = "pathname already exists (not necessarily as a directory). This includes the case where pathname is a symbolic link, dangling or not.";
				break;
			case EFAULT:
				error_message = "pathname points outside your accessible address space.";
				break;
			case ELOOP:
				error_message = "Too many symbolic links were encountered in resolving pathname.";
				break;
			case EMLINK:
				error_message = "The number of links to the parent directory would exceed LINK_MAX.";
				break;
			case ENAMETOOLONG:
				error_message = "pathname was too long.";
				break;
			case ENOENT:
				error_message = "A directory component in pathname does not exist or is a dangling symbolic link.";
				break;
			case ENOMEM:
				error_message = "Insufficient kernel memory was available.";
				break;
			case ENOSPC:
				error_message = "The device containing pathname has no room for the new directory or the user's disk quota is exhausted.";
				break;
			case ENOTDIR:
				error_message = "A component used as a directory in pathname is not, in fact, a directory.";
				break;
			case EPERM:
				error_message = "The filesystem containing pathname does not support the creation of directories.";
				break;
			case EROFS:
				error_message = "pathname refers to a file on a read-only filesystem.";
				break;
			default:
				error_message = "Unknown error occurred while creating directory.";
				break;
			}
			throw ofstream::failure(error_message);
		}
	}
};

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_UTILS_MAKEDIRECTORY_H_ */
