/*
 * PotentialGridDX.cpp
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#include "../exceptions/ParsingOpenDXException.h"
#include "PotentialGridDX.h"
#include <fstream>
#include <math.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

/**
 * Constructor for the potentials grid derived from the APBS package.
 * @param openDX_filename		name of the file containing the electrostatic
 * 								potentials calculated with the APBS tool
 */
PotentialGridDX::PotentialGridDX(string const & openDX_filename) {

	Vcompare = pow(10,-1*(VGRID_DIGITS - 2));

	std::ifstream file_stream;
	file_stream.open(openDX_filename);
	if (!file_stream.is_open()) {
		struct stat st = { 0 };
		stat(openDX_filename.c_str(), &st);
		string error_message;
		switch (errno) {
		case EACCES:
			error_message = "Search permission is denied for one of the directories in the path prefix of path.";
			break;
		case EBADF:
			error_message = "File descriptor is bad.";
			break;
		case EFAULT:
			error_message = "Pathname points outside your accessible address space.";
			break;
		case ELOOP:
			error_message = "Too many symbolic links were encountered in resolving pathname.";
			break;
		case ENAMETOOLONG:
			error_message = "Path is too long.";
			break;
		case ENOENT:
			error_message = "The file or a directory component in pathname does not exist or is a dangling symbolic link.";
			break;
		case ENOMEM:
			error_message = "Insufficient kernel memory was available.";
			break;
		case ENOTDIR:
			error_message = "A component used as a directory in pathname is not, in fact, a directory.";
			break;
		case EOVERFLOW:
			error_message = "Path refers to a file whose size, inode number, or number of blocks cannot be represented in, respectively, the types off_t, ino_t, or blkcnt_t. This error can occur when, for example, an application compiled on a 32-bit platform without -D_FILE_OFFSET_BITS=64 calls stat() on a file whose size exceeds (1<<31)-1 bytes.";
			break;
		default:
			error_message = "Unknown error occurred while opening file.";
			break;
		}
		throw std::ifstream::failure("Failed opening input DX file: " + openDX_filename + " " + error_message);
	}
	int numscan;
	int tx, ty, tz, n;
	double td1, td2, td3;
	char line[255];
	while (file_stream.getline(line, 255)) {
		/* if blank line or comment continue */
		if (line[0] == '#' || line[0] == '\n' || line[0] == '\0')
			continue;
		else
			break;
	}
	numscan = sscanf(line, "object 1 class gridpositions counts %d %d %d",
			&nx, &ny, &nz);
	if (numscan == 3)
		file_stream.getline(line, 255);
	else {
		std::stringstream error;
		error << "Malformed input OpenDX file: " << openDX_filename << ".";
		std::stringstream cause;
		cause << "Number of parsed fields is: " << numscan
				<< ". Should be 3.";
		throw ParsingOpenDXException(error.str(), "PotentialGridDX()",
				"Incorrect OpenDX file format. " + cause.str());
	}
	numscan = sscanf(line, "origin %lf %lf %lf", &xmin, &ymin, &zmin);
	if (numscan == 3)
		file_stream.getline(line, 255);
	else {
		std::stringstream error;
		error << "Malformed input OpenDX file: " << openDX_filename << ".";
		std::stringstream cause;
		cause << "Number of parsed fields is: " << numscan
				<< ". Should be 3.";
		throw ParsingOpenDXException(error.str(), "PotentialGridDX()",
				"Incorrect OpenDX file format. " + cause.str());
	}

	numscan = sscanf(line, "delta %lf %*s %*s", &hx);
	if (numscan == 1)
		file_stream.getline(line, 255);
	else {
		std::stringstream error;
		error << "Malformed input OpenDX file: " << openDX_filename << ".";
		std::stringstream cause;
		cause << "Number of parsed fields is: " << numscan
				<< ". Should be 1.";
		throw ParsingOpenDXException(error.str(), "PotentialGridDX()",
				"Incorrect OpenDX file format. " + cause.str());
	}
	numscan = sscanf(line, "delta %*s %lf %*s", &hy);
	if (numscan == 1)
		file_stream.getline(line, 255);
	else {
		std::stringstream error;
		error << "Malformed input OpenDX file: " << openDX_filename << ".";
		std::stringstream cause;
		cause << "Number of parsed fields is: " << numscan
				<< ". Should be 1.";
		throw ParsingOpenDXException(error.str(), "PotentialGridDX()",
				"Incorrect OpenDX file format. " + cause.str());
	}
	numscan = sscanf(line, "delta %*s %*s %lf", &hz);
	if (numscan == 1)
		file_stream.getline(line, 255);
	else {
		std::stringstream error;
		error << "Malformed input OpenDX file: " << openDX_filename << ".";
		std::stringstream cause;
		cause << "Number of parsed fields is: " << numscan
				<< ". Should be 1.";
		throw ParsingOpenDXException(error.str(), "PotentialGridDX()",
				"Incorrect OpenDX file format. " + cause.str());
	}
	numscan = sscanf(line, "object 2 class gridconnections counts %d %d %d", &tx, &ty, &tz);
	if (numscan == 3 && tx == nx && ty == ny && tz == nz)
		file_stream.getline(line, 255);
	else if (numscan != 3) {
		std::stringstream error;
		error << "Malformed input OpenDX file: " << openDX_filename << ".";
		std::stringstream cause;
		cause << "Number of parsed fields is: " << numscan
				<< ". Should be 3.";
		throw ParsingOpenDXException(error.str(), "PotentialGridDX()",
				"Incorrect OpenDX file format. " + cause.str());
	}
	else {
		std::stringstream error;
		error << "Malformed input OpenDX file: " << openDX_filename << ".";
		std::stringstream cause;
		cause << "object 1 class gridpositions counts and object 2 class gridconnections counts do not match.";
		throw ParsingOpenDXException(error.str(), "PotentialGridDX()",
				"Incorrect OpenDX file format. " + cause.str());
	}
	numscan = sscanf(line, "object 3 class array type double rank 0 items %d data follows", &n);
	if (numscan == 1 && n == nx * ny * nz)
		file_stream.getline(line, 255);
	else if (numscan != 1) {
		std::stringstream error;
		error << "Malformed input OpenDX file: " << openDX_filename << ".";
		std::stringstream cause;
		cause << "Number of parsed fields is: " << numscan
				<< ". Should be 1.";
		throw ParsingOpenDXException(error.str(), "PotentialGridDX()",
				"Incorrect OpenDX file format. " + cause.str());
	}
	else {
		std::stringstream error;
		error << "Malformed input OpenDX file: " << openDX_filename << ".";
		std::stringstream cause;
		cause << "The number of object 3 class array type double rank 0 items does not match the parsed dimensions of the grid.";
		throw ParsingOpenDXException(error.str(), "PotentialGridDX()",
				"Incorrect OpenDX file format. " + cause.str());
	}
    xmax = xmin + (nx-1)*hx;
    ymax = ymin + (ny-1)*hy;
    zmax = zmin + (nz-1)*hz;

	data = new double[n];

	int counter = 0;
	while (counter < n) {
		numscan = sscanf(line, "%lf %lf %lf", &td1, &td2, &td3);
		if (numscan == 3) {
			data[counter] = td1;
			data[counter + 1] = td2;
			data[counter + 2] = td3;
		} else if (numscan == 2) {
			data[counter] = td1;
			data[counter + 1] = td2;
		} else if (numscan == 1) {
			data[counter] = td1;
		} else {
			std::stringstream error;
			error << "Malformed input OpenDX file: " << openDX_filename << ".";
			std::stringstream cause;
			cause << "The data ended unexpectedly. " << n
					<< " entries expected, but there were only " << counter
					<< ".";
			throw ParsingOpenDXException(error.str(), "PotentialGridDX()",
					"Incorrect OpenDX file format. " + cause.str());		}
		counter += numscan;
		file_stream.getline(line, 255);
	}

	file_stream.close();
}
/**
 * Destructor of the class.
 */
PotentialGridDX::~PotentialGridDX() {
	delete[] data;
	data = NULL;
}
