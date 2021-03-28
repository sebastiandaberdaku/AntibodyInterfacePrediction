/*
 * CommandLineParser.cpp
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

#include "../utils/disclaimer.h"
#include "CommandLineParser.h"
#include "CustomValidators.h"
/**
 * This method parses the input command line options.
 * @param argc		input argc reference, i.e. the first argument of the main method
 * @param argv		input argv reference, i.e. the second argument of the main method
 *
 * @return 			true if all arguments are parsed correctly, false otherwise.
 */
bool CommandLineParser::parseCommandLine(int const & argc,
		char const * const argv[], float & patchRadius, float & minCenterDist_noninterface,
		float & minCenterDist_interface,
		float & probeRadius, float & resolution, string & inname1,
		string & inname2, string & outname1,
		string & outname2, string & inname_radii, bool & no_hydrogen,
		bool & no_hetatm, int & maxOrder, float & interface_distance,
		bool & help, bool & version, bool & surf_description, bool & license) {

	//------------------------------- command line options --------------------------------------
	options_description description("Paratope prediction - Training samples generator. Usage");
	description.add_options()("help,h", "Display this brief help message.")
	/*
	 * This token is used to specify the name of the first input file. It is mandatory.
	 * The receptor command line token has no option name. The command line tokens which have
	 * no option name are called "positional options" in the boost::program_options library.
	 */
	("receptor", value<input_PDB_PQR_filename>()->required(), "Name of the Antibody's PDB file.")
	/*
	 * This token is used to specify the name of the second input file. It is mandatory.
	 * The ligand command line token has no option name. The command line tokens which have
	 * no option name are called "positional options" in the boost::program_options library.
	 */
	("ligand", value<input_PDB_PQR_filename>()->required(), "Name of the Antigen's PDB file.")
	/*
	 * This token is used to specify the input atom radii filename.
	 */
	("atom_radii", value<filename>(), "File containing the radius information of each atom. If not specified, the default CHARMM22 radius values will be used for PDB files. It is ignored if specified for PQR files.")
	/*
	 * The --no_hetatm flag is used to include HETATM records in the surface computation.
	 */
	("no_hetatm", "Ignore HETATM records in the surface computation.")
	/*
	 * The --no_hydrogen flag is used to include hydrogen atoms in the surface computation.
	 */
	("no_hydrogen", "Ignore hydrogen atoms in the surface computation.\nNOTE: X-ray crystallography cannot resolve hydrogen atoms \
in most protein crystals, so in most PDB files, hydrogen atoms are absent. Sometimes hydrogens are added by modeling. \
Hydrogens are always present in PDB files resulting from NMR analysis, and are usually present in theoretical models.")
	/*
	 * The -p (--probe_radius) flag is used to specify the probe radius. Default is 1.4. Optional.
	 */
	("probe_radius,p", value<probe_radius>()->default_value(1.4), "Probe radius (in Å), floating point number in (0, 2.0] (default is 1.4Å).")
	/*
	 * The -R(--patch_radius) flag is used to specify the radius of the sphere which determines the surface patch.
	 * Default is 6. Optional.
	 */
	("patch_radius,R", value<patch_rd>()->default_value(6.0), "Patch radius (in Å), positive floating point (default is 6.0Å).")
	/*
	 * The -i(--interface_distance) flag is used to specify the threshold for the interface computation.
	 * Default is 4.5. Optional.
	 */
	("interface_distance,i", value<patch_rd>()->default_value(6.0), "Threshold distance for the interface determination (in Å), positive floating point (default is 6.0Å).")
	/*
	 * The -d (--patch_dist) flag is used to specify the minimum distance between patch centers.
	 * Default is 1.0. Optional.
	 */
	("patch_dist_i", value<patch_rd>()->default_value(1.0), "Minimum distance between interface patch centers (in Å), positive floating point (default is 1.0Å).")
	("patch_dist_n", value<patch_rd>()->default_value(3.5), "Minimum distance between non-interface patch centers (in Å), positive floating point (default is 3.5Å).")

	/*
	 * The -r (--resolution) flag is used to specify the resolution of the voxel grid. Default is 4. Optional.
	 */
	("resolution,r", value<resolution_param>()->default_value(4.0), "Resolution factor, positive floating point (default is 4.0). This value's cube determines the number of voxels per Å³.")
	/*
	 * The -N (--max_order) flag is used to specify the maximum order for the Zernike descriptor of the surface patches.
	 * Default is 20. Optional.
	 */
	("max_order,N", value<max_order>()->default_value(20), "Maximum Zernike descriptor order (default is 20)")
	/*
	 * The -l (--license) flag is used to view the program's license
	 * information.
	 */
	("license,l", "View license information.")
	/*
	 * The -v (--version) flag is used to view the program's version
	 * information.
	 */
	("version,v", "Display the version number");
	/*
	 * The receptorPDB and ligandPDB options must be declared as positional.
	 */
	positional_options_description p;
	p.add("receptor", 1);
	p.add("ligand", 1);

	variables_map vm;

	try {
		//--------------------------------parsing command line options------------------------------
		/*
		 * And it is finally specified when parsing the command line.
		 */
		store(command_line_parser(argc, argv).options(description).positional(p).run(), vm);

		if (vm.count("help")) {
			cout << description;
			help = true;
			return true;
		}
		if (vm.count("version")) {
			cout << "Program: " << argv[0] << ", version: " << PROGRAM_VERSION << "\n";
			version = true;
			return true;
		}
		if (vm.count("license")) {
			DISCLAIMER
			license = true;
			return true;
		}
		/*
		 * notify throws exceptions so we call it after the above checks
		 */
		notify(vm);

		/* initializing variables */

		inname1 = vm["receptor"].as<input_PDB_PQR_filename>().filename;
		int lastindex = inname1.find_last_of(".");
		outname1 = inname1.substr(0, lastindex);

		inname2 = vm["ligand"].as<input_PDB_PQR_filename>().filename;
		lastindex = inname2.find_last_of(".");
		outname2 = inname2.substr(0, lastindex);

		probeRadius = vm["probe_radius"].as<probe_radius>().p;
		patchRadius = vm["patch_radius"].as<patch_rd>().p;

		minCenterDist_interface = vm["patch_dist_i"].as<patch_rd>().p;
		minCenterDist_noninterface = vm["patch_dist_n"].as<patch_rd>().p;

		interface_distance = vm["interface_distance"].as<patch_rd>().p;
		resolution = vm["resolution"].as<resolution_param>().r;
		maxOrder = vm["max_order"].as<max_order>().n;

		if (vm.count("atom_radii")) {
			inname_radii = vm["atom_radii"].as<filename>().fname;
		} else {
			inname_radii = "CHARMM22";
		}
		no_hetatm = vm.count("no_hetatm");
		no_hydrogen = vm.count("no_hydrogen");
	} catch (error const & e) {
		cerr << "error: " << e.what() << "\n";
		cerr << description << "\n";
		return false;
	}
	return true;
};


