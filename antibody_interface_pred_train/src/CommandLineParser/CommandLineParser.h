/*
 * CommandLineParser.h
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

#ifndef COMMANDLINEPARSER_H_
#define COMMANDLINEPARSER_H_

#include <string>

using namespace std;

/**
 * This class provides a simple interface for accessing command line arguments.
 */
class CommandLineParser {
public:
	static bool parseCommandLine(int const & argc, char const * const argv[],
			float & patchRadius, float & minCenterDist_noninterface,
			float & minCenterDist_interface, float & probeRadius,
			float & resolution, string & inname1, string & inname2,
			string & outname1, string & outname2,
			string & inname_radii, bool & no_hydrogen,
			bool & no_hetatm, int & maxOrder, float & interface_distance,
			bool & help, bool & version,
			bool & surf_description, bool & license);
};

#endif /* COMMANDLINEPARSER_H_ */
