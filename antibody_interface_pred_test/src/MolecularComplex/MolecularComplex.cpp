#include "MolecularComplex.h"
#include "../Geometry/normalizeGrid.h"
#include "../utils/makeDirectory.h"
/*
 * MolecularComplex.cpp
 *
 *  Created on: 03/nov/2014
 *      Author: sebastian
 */

MolecularComplex::MolecularComplex(Molecule const & receptor, Molecule const & ligand, float interface_distance /** 6.0Å distance used to determine if an atom belongs to the interface*/) :
	receptor(&receptor), ligand(&ligand), dockingPose(NULL) {

	float s_searchRad = pow(interface_distance, 2.0);

	for (auto const & atm : *receptor.atoms) {

		vector<pair<size_t, float> > ret_matches;
		size_t n = ligand.atoms_tree.radiusSearch(atm, s_searchRad, ret_matches);
		if (n > 0) {
			interface_receptor_residues.insert(residue(atm.chain_ID, atm.atom_number));
		}
	}

	receptor_interface = array3D(receptor.length, array2D(receptor.width, array1D(receptor.height, 0)));

	for (int i  = 0; i < receptor.length; ++i) {
		for (int j = 0; j < receptor.width; ++j) {
			for (int k = 0; k < receptor.height; ++k) {
				if (receptor.surface->surface->getVoxel(i, j, k)) {
					point3D cp(i, j, k);
					cp /= receptor.surface->resolution;
					cp -= receptor.translation;

					float dist;
					size_t n_atm_id = receptor.atoms_tree.nnSearch(cp, dist);
					residue res((*receptor.atoms)[n_atm_id].chain_ID, (*receptor.atoms)[n_atm_id].atom_number);
					if (interface_receptor_residues.count(res) > 0) {
						receptor_interface[i][j][k] = 1;
					}
				}
			}
		}
	}


};
MolecularComplex::MolecularComplex(Molecule const & receptor, Molecule const & ligand, DockingPose const & p) :
	receptor(&receptor), ligand(&ligand), dockingPose(&p) {
};

MolecularComplex::~MolecularComplex() {
};
/**bool MolecularComplex::calculateScores(double & hydrophobicity, double & bindingEnergy,
		double & CoulombicPotential, double & ACE) {
	vector<atom> ligand_atoms(*(ligand->atoms));
	transformAtoms(ligand_atoms, *dockingPose);
	if (compenetrationAllowed(ligand_atoms)) {
		hydrophobicity = getHydrophobicEnergy(ligand_atoms);
		bindingEnergy = getBindingEnergy(ligand_atoms);
		CoulombicPotential = getCoulombicEnergy(ligand_atoms);
		ACE = getAtomicContactEnergy(ligand_atoms);
		return true;
	}
	hydrophobicity = min_double;
	bindingEnergy = min_double;
	CoulombicPotential = min_double;
	ACE = min_double;
	return false;
};*/
/**
 * Discard poses which yield high intermolecular penetration.
 *
 * By allowing some intermolecular penetrations we implicitly take into
 * account a certain extent of conformational flexibility. The only
 * solutions which are discarded are those in which ligand atoms fall
 * into the "core" of the receptor protein. "Core" atoms are those
 * atoms that have zero accessible surface area. Solutions with ligand
 * atom centers which invade the outer shell of the molecular
 * representation are retained.
 * @return true if the compenetration is acceptable, false otherwise
 */
//bool MolecularComplex::compenetrationAllowed() {
//		vector<atom> ligand_atoms(ligand->pqrModel->atomsInModel);
//		transformAtoms(ligand_atoms, *dockingPose);
//		return compenetrationAllowed(ligand_atoms);
//};
/** Prints the 3D voxelized representation to file using the PCD
 * (Point Cloud Data) file format.
 *
 * Each PCD file contains a header that identifies and declares
 * certain properties of the point cloud data stored in the file.
 * The header of a PCD must be encoded in ASCII.
 * Storing point cloud data in both a simple ascii form with each
 * point on a line, space or tab separated, without any other
 * characters on it, as well as in a binary dump format, allows
 * us to have the best of both worlds: simplicity and speed,
 * depending on the underlying application. The ascii format
 * allows users to open up point cloud files and plot them using
 * standard software tools like gnuplot or manipulate them using
 * tools like sed, awk, etc.
 *
 * For a detailed description of the PCD (Point Cloud Data) file
 * format specification see:
 * http://pointclouds.org/documentation/tutorials/pcd_file_format.php
 *
 * \param filename	Name of the output file. The '.pcd' extension
 * 					is added automatically.
 * \param algorithmType
 * \throws ofstream::failure
 */
/**void MolecularComplex::outputSurfacePCDModel(string const & filename) {
	time_t time = clock();
	uint16_t length, width, height;
	point3D translation;
	vector<atom> complex_atoms(ligand->pqrModel->atomsInModel);
	transformAtoms(complex_atoms, *dockingPose);
	complex_atoms.insert(complex_atoms.begin(), receptor->pqrModel->atomsInModel.begin(), receptor->pqrModel->atomsInModel.end());

	cout << "Initializing parameters for the merged PDB model.\n";
	MolecularSurface::boundingBox(complex_atoms, receptor->surface->probeRadius,
			receptor->surface->resolution, length, width, height, translation, Molecule::max_atm_radius, Molecule::min_atm_radius);
	MolecularSurface surface(complex_atoms, receptor->surface->probeRadius, receptor->surface->resolution,
			length, width, height, translation);

	cout << "Creating the space-filling model\n";
	surface.createCPKModel();
	cout << "Filling internal cavities\n";
	surface.fillInternalCavities();
	surface.buildSurface();
	cout << "Calculating the Euclidean Distance Transform using \n";
	cout << "speed-optimized data structures\n";
	surface.fastRegionGrowingEDT();
	cout << "Surface calculation time:\t" << (double)(clock() - time) / CLOCKS_PER_SEC << " seconds.\n";
	// Output PCD model to file 
	cout << "Output surface PCD model\n";
	surface.outputSurfacePCDModel(filename);
};*/
void MolecularComplex::outputReceptorInterfacePCDModel(string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + filename + ".pcd");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream << std::setprecision(11);
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < receptor->length; ++i) {
		for (int j = 0; j < receptor->width; ++j) {
			for (int k = 0; k < receptor->height; ++k) {
				if (receptor_interface[i][j][k] > 0)
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}


	union {
		uint8_t rgb_byte[4];
		float rgb_flt;
	} color_RGB;

	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n" << "FIELDS x y z rgb\n" << "SIZE 4 4 4 "<< sizeof(color_RGB) <<"\n"
			<< "TYPE F F F F\n" << "COUNT 1 1 1 1\n" << "WIDTH "
			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
			<< "\n" << "DATA ascii";
	float pDX; /**< current potential */
	int x, y, z;
	uint8_t r, g, b;

	while (!surfaceVoxels.empty()) {
		x = surfaceVoxels.front().ix;
		y = surfaceVoxels.front().iy;
		z = surfaceVoxels.front().iz;

		pDX = 1.0; //r_positiveDXPotentials[x][y][z] - r_negativeDXPotentials[x][y][z];
		if (pDX <= 0) {
			r = 255;
			g = (uint8_t) floor((1.0 + pDX) * 255);
			b = (uint8_t) floor((1.0 + pDX) * 255);
		} else {
			r = (uint8_t) floor((1.0 - pDX) * 255);
			g = (uint8_t) floor((1.0 - pDX) * 255);
			b = 255;
		}

		color_RGB.rgb_byte[0] = b;
		color_RGB.rgb_byte[1] = g;
		color_RGB.rgb_byte[2] = r;
		color_RGB.rgb_byte[3] = 0;

		file_stream << "\n" << x / receptor->surface->resolution - receptor->translation.x << " "
				<< y / receptor->surface->resolution - receptor->translation.y << " "
				<< z / receptor->surface->resolution - receptor->translation.z << " "
				<< color_RGB.rgb_flt;
		surfaceVoxels.pop();
	}
	file_stream.close();
}

//void MolecularComplex::outputLigandInterfacePCDModel(string const & filename) {
//	std::ofstream file_stream;
//	makeDirectory("./output");
//	file_stream.open("./output/" + filename + ".pcd");
//	if (!file_stream.is_open()) {
//		throw std::ofstream::failure("Error opening output file.");
//	}
//	file_stream << std::setprecision(11);
//	std::queue<voxel> surfaceVoxels;
//	for (int i = 0; i < ligand->length; ++i) {
//		for (int j = 0; j < ligand->width; ++j) {
//			for (int k = 0; k < ligand->height; ++k) {
//				if (ligand_interface[i][j][k])
//					surfaceVoxels.push(voxel(i, j, k));
//			}
//		}
//	}
//
//	union {
//		uint8_t rgb_byte[4];
//		float rgb_flt;
//	} color_RGB;
//
//	/* File header */
//	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
//			<< "VERSION .7\n" << "FIELDS x y z rgb\n" << "SIZE 4 4 4 "<< sizeof(color_RGB) <<"\n"
//			<< "TYPE F F F F\n" << "COUNT 1 1 1 1\n" << "WIDTH "
//			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
//			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
//			<< "\n" << "DATA ascii";
//	float pDX; /**< current potential */
//	int x, y, z;
//	uint8_t r, g, b;
//
//	while (!surfaceVoxels.empty()) {
//		x = surfaceVoxels.front().ix;
//		y = surfaceVoxels.front().iy;
//		z = surfaceVoxels.front().iz;
//
//		pDX = l_positiveDXPotentials[x][y][z] - l_negativeDXPotentials[x][y][z];
//		if (pDX <= 0) {
//			r = 255;
//			g = (uint8_t) floor((1.0 + pDX) * 255);
//			b = (uint8_t) floor((1.0 + pDX) * 255);
//		} else {
//			r = (uint8_t) floor((1.0 - pDX) * 255);
//			g = (uint8_t) floor((1.0 - pDX) * 255);
//			b = 255;
//		}
//
//		color_RGB.rgb_byte[0] = b;
//		color_RGB.rgb_byte[1] = g;
//		color_RGB.rgb_byte[2] = r;
//		color_RGB.rgb_byte[3] = 0;
//
//		file_stream << "\n" << x / ligand->surface->resolution - ligand->translation.x << " "
//				<< y / ligand->surface->resolution - ligand->translation.y << " "
//				<< z / ligand->surface->resolution - ligand->translation.z << " "
//				<< color_RGB.rgb_flt;
//		surfaceVoxels.pop();
//	}
//	file_stream.close();
//
//}
