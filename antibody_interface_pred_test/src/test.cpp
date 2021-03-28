/*
 * validation.cpp
 *
 *  Created on: 6/mag/2016
 *      Author: sebastian
 */

#include "CommandLineParser/CommandLineParser.h"
#include "DockingMethods/DockingMethods.h"
#include "exceptions/ParsingOpenDXException.h"
#include "exceptions/ParsingPQRException.h"
#include "hydrophobicity/hydrophobicity.h"
#include "MolecularComplex/MolecularComplex.h"
#include "Molecule/Molecule.h"
#include "utils/disclaimer.h"
#include "utils/elapsedTime.h"
#include "utils/makeDirectory.h"
#include "Zernike/BoundingSphere.h"
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <utility>

#include <functional>

#include "test.h"

using namespace std;

int main (int argc, char* argv[]) {
//---------------------------variables and parameters------------------------------------
	float patchRadius; // patch sphere radius
	float minCenterDist; // minimum distance between patch centers
	float probeRadius; // probe sphere radius
	float resolution; // resolution^3 = #voxels/Å^3

	float interface_distance;

	string inname_receptor; // input filename
	string inname_ligand; // input filename

	string outname_receptor; // output filename
	string outname_ligand; // output filename

	string inname_radii;

	bool no_hydrogen;
	bool no_hetatm;

	bool help = false; // if true, print the help message
	bool version = false; // if true, print the program version
	bool surf_description = false; // if true, print the three surface descriptions
	bool license = false; // if true, print the license information


	int maxOrder; // max Zernike Descriptor order - N in paper

    auto const t_start = chrono::high_resolution_clock::now();

	if (!CommandLineParser::parseCommandLine(argc, argv, patchRadius,
			minCenterDist, probeRadius, resolution, inname_receptor,
			inname_ligand, outname_receptor, outname_ligand,
			inname_radii, no_hydrogen, no_hetatm, maxOrder,
			interface_distance, help, version, license))
		return EXIT_FAILURE;
	if (help || version || surf_description || license)
		return EXIT_SUCCESS;
//-----------------------------print config---------------------------------------------
	PROGRAM_INFO
	/* summary of the parsed parameters */
	cout << "The specification is: \n" << "input filenames: " << inname_receptor << " and " << inname_ligand << "\n";
	cout << "SES computation algorithm:\tRegion Growing EDT with speed-optimized \n";
	cout << "\t\tdata structures, \n";
	cout << "probe radius:\t" << probeRadius << "Å, \n";
	cout << "resolution:\t" << pow((double) resolution, 3.0)
			<< " voxels per Å³, \n";
	cout << "interface distance threshold: " << interface_distance << "Å, \n";

	if (no_hetatm)
		cout << "include HETATM records: no\n";
	else
		cout << "include HETATM records: yes\n";
	if (no_hydrogen)
		cout << "include hydrogen atoms: no\n";
	else
		cout << "include hydrogen atoms: yes\n";

	string extension_receptor = inname_receptor.substr(inname_receptor.find_last_of(".") + 1);
	string extension_ligand = inname_ligand.substr(inname_ligand.find_last_of(".") + 1);
	if(extension_ligand == "pdb" || extension_ligand == "PDB"
			|| extension_receptor == "pdb" || extension_receptor == "PDB") {
		cout << "atomic radii: " << inname_radii << "\n";
	}
	cout << "patch radius:\t" << patchRadius << "Å, \n";
	cout << "minimum distance between patch centers:\t" << minCenterDist
			<< "Å, \n";
	cout << "maximum Zernike descriptor order:\t" << maxOrder << ". \n";
	cout << "**************************************************\n";
//-----------------------------computation--------------------------------------------
	try {
		Molecule *receptor, *ligand;
#pragma omp parallel sections
		{
#pragma omp section
			{
				receptor = new Molecule(patchRadius, minCenterDist, probeRadius,
						resolution, inname_receptor,
						outname_receptor, inname_radii, maxOrder, no_hydrogen,
						no_hetatm);
			}
#pragma omp section
			{
				ligand = new Molecule(patchRadius, minCenterDist, probeRadius,
						resolution, inname_ligand,
						outname_ligand, inname_radii, maxOrder, no_hydrogen,
						no_hetatm, false);
			}
		}
		MolecularComplex complex(*receptor, *ligand, interface_distance);

		vector<CompactPatchDescriptor> receptor_descriptors;

		receptor->calculateDescriptors(inname_receptor, minCenterDist, patchRadius, maxOrder, complex.receptor_interface, receptor_descriptors);

		export_descriptors(receptor_descriptors, maxOrder, inname_receptor.substr(0, inname_receptor.length() - 4));


		cout << "**************************************************\n";
		cout << "Total calculation time:\t" << elapsedTime(t_start, chrono::high_resolution_clock::now()) << "\n";
		cout << "**************************************************\n";
		delete receptor;
		delete ligand;
	} catch (ParsingPQRException const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (ParsingOpenDXException const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (fstream::failure const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (out_of_range const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (invalid_argument const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (logic_error const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (runtime_error const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
