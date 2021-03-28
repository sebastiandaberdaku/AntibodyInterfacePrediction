/*
 * Molecule.cpp
 *
 *  Created on: 07/set/2014
 *      Author: sebastian
 */

#include "../ASP/AtomicSolvationParameters.h"
#include "../hydrophobicity/hydrophobicity.h"
#include "../utils/doubleCompare.h"
#include "../utils/makeDirectory.h"
#include "Molecule.h"
#include <omp.h>
Molecule::Molecule(float patchRadius, float minCenterDist,
		float probeRadius, float resolution, string const & inname,
		string const & outname, string const & inname_radii,
		int maxOrder, bool no_hydrogen, bool no_hetatm,  bool receptor) :
			ASA(0), dG(0), SEV(0), pqrModel(NULL), pdbModel(NULL), surface(NULL) {

	max_atm_radius = min_float;
	min_atm_radius = max_float;
	n_sphere_points = 96;

	string extension = inname.substr(inname.find_last_of(".") + 1);
	if(extension == "pqr" || extension == "PQR") {
		cout << "Process " << omp_get_thread_num() << ": " << "Loading PQR file: "<< inname <<"\n";
		pqrModel = new PQRModel(inname, no_hydrogen, no_hetatm);
		atoms = &pqrModel->atomsInModel;
	} else if (extension == "pdb" || extension == "PDB") {
		cout << "Process " << omp_get_thread_num() << ": " << "Loading PDB file: "<< inname <<"\n";
		pdbModel = new PDBModel(inname, inname_radii, no_hydrogen, no_hetatm);
		atoms = &pdbModel->atomsInModel;
	}

	atoms_tree = kdtree_atom(atoms);

	if (receptor) {

		cout << "Process " << omp_get_thread_num() << ": "
				<< "Initializing parameters for " << inname << ".\n";
		cout << "Process " << omp_get_thread_num() << ": "
				<< "Number of atoms in model: " << atoms->size() << "\n";
		MolecularSurface::boundingBox(*atoms, probeRadius, patchRadius,
				resolution, length, width, height, translation, max_atm_radius,
				min_atm_radius);
		auto t_start = chrono::high_resolution_clock::now();

		surface = new MolecularSurface(*atoms, probeRadius, resolution, length,
				width, height, translation);
		cout << "Process " << omp_get_thread_num() << ": "
				<< "Creating the space-filling model\n";
		surface->createCPKModel();
		cout << "Process " << omp_get_thread_num() << ": "
				<< "Filling internal cavities\n";
		surface->fillInternalCavities();
		cout << "Process " << omp_get_thread_num() << ": "
				<< "Building SA surface\n";
		surface->buildSurface_old();

//	cout << "Process " << omp_get_thread_num() << ": " << "Calculating Accessible Surface Area\n";
//	ASA = surface->calculateAccessibleSurfaceArea(perAtomASA, per_atom_SA_points, n_sphere_points);
//	cout << "Process " << omp_get_thread_num() << ": " << "Calculating Solvation Energy\n";
//	dG = surface->calculateSolvationEnergy(AtomicSolvationParameters::asp_zhou_2, perAtomASA);

//	size_t num_atoms = atoms->size();
//	for (size_t ii = 0; ii < num_atoms; ++ii) {
//		if (doubleCompare(perAtomASA[ii], 0))
//			core_atoms.push_back((*atoms)[ii]);
//		else
//			outer_atoms.push_back((*atoms)[ii]);
//	}
//	core_atoms_tree = kdtree_atom(&core_atoms);
//	outer_atoms_tree = kdtree_atom(&outer_atoms);

		cout << "Process " << omp_get_thread_num() << ": "
				<< "Calculating the Euclidean Distance Transform using \n";
		cout << "Process " << omp_get_thread_num() << ": "
				<< "speed-optimized data structures\n";
		surface->fastRegionGrowingEDT();
		/*
		 cout << "Process " << omp_get_thread_num() << ": " << "Calculating the Solvent-Excluded Volume\n";
		 SEV = surface->calculateSolventExcludedVolume();
		 */
//	cout << "Process " << omp_get_thread_num() << ": " << "Mapping the APBS-calculated potentials on the \nmolecular surface.\n";
//	surface->calculateDXPotentials(openDX);
//	cout << "Process " << omp_get_thread_num() << ": " << "Calculating the hydrophobicity map.\n";
//	surface->calculateHydrophobicity(hphob_kyte, 15);
//	surface->calculateAtomHydrophobicity(hphob_kapcha);
		cout << "Process " << omp_get_thread_num() << ": "
				<< "Calculating HQI8 indices.\n";
		surface->calculateHQI8();

		auto t_end = chrono::high_resolution_clock::now();
		int t_ms =
				chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
		cout << "Process " << omp_get_thread_num() << ": "
				<< "Surface calculation time:\t" << t_ms / 1000.0
				<< " seconds.\n";

		/* Output PCD model to file */
//	cout << "Process " << omp_get_thread_num() << ": " << "Output surface PCD model\n";
//	surface->outputSurfacePCDModel(outname);
//	surface->outputDXPotentials(outname);
//	cout << "Process " << omp_get_thread_num() << ": " << "Output hydrophobicity PCD model\n";
//	surface->outputHydrophobicity(outname);

//	t_start = chrono::high_resolution_clock::now();
//	cout << "Process " << omp_get_thread_num() << ": " << "Extracting patch centers.\n";
//	surface->extractPatchCenters(minCenterDist);
//	surface->outputPatchCentersPCDModel(outname);
//	cout << "Process " << omp_get_thread_num() << ": " << "Calculating Surface Descriptors\n";
//	surface->calculateSurfaceDescriptors(patchRadius, max_atm_radius, maxOrder, descriptors);
//	CompactPatchDescriptor::set_curvature_type(descriptors);
//	t_end = chrono::high_resolution_clock::now();
//	t_ms = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
//	outputMoleculeDescriptors(outname);
//	cout << "Process " << omp_get_thread_num() << ": " << "Patch extraction and descriptor \ncalculation time:\t" << t_ms / 1000.0 << " seconds.\n"<<endl;
	}
}

Molecule::~Molecule() {
	if (surface != NULL) {
		delete surface;
		surface = NULL;
	}
	if (pqrModel != NULL) {
		delete pqrModel;
		pqrModel = NULL;
	}
	if (pdbModel != NULL) {
		delete pdbModel;
		pdbModel = NULL;
	}
}

void Molecule::outputMoleculePQR(string const & filename) {
	this->pqrModel->outputPQRFile(filename);
}

void Molecule::outputMoleculePQR(string const & filename, DockingPose const & p) {
	PQRModel model(*this->pqrModel);
	model.transform(p);
	model.outputPQRFile(filename);
}
void Molecule::outputMoleculeDescriptors(string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	makeDirectory("./output/descriptors");
	file_stream.open("./output/descriptors/" + filename + "-descriptors.txt");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	/* File header */
	file_stream << descriptors.size();//"Zernike descriptors for: " << filename;
	for (auto const & cpd : descriptors) {
		file_stream << "\n" << cpd;
	}
	file_stream.close();
}

