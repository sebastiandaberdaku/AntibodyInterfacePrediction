/*
 * Molecule.h
 *
 *  Created on: 07/set/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_MOLECULE_MOLECULE_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_MOLECULE_MOLECULE_H_

#include "../DockingMethods/DockingPose.h"
#include "../kdtree/kdtree.h"
#include "../MolecularSurface/MolecularSurface.h"
#include "../PQR/PQRModel.h"
#include "../PDB/PDBModel.h"
#include <string.h>
#include <chrono>
#include <ctime>
#include <cstdlib>      // std::rand, std::srand

using namespace std;

class Molecule {
public:
	MolecularSurface * surface;
	PQRModel * pqrModel;
	PDBModel * pdbModel;
	uint16_t length, width, height; // bounding box dimensions
	point3D translation; // translation vector
	double SEV; /**< Solvent-Excluded Volume */
	double ASA; /**< Accessible Surface Area */
	double dG; /**< Solvation Energy */
	vector<double> perAtomASA;
	vector<vector<point3D>> per_atom_SA_points;

	vector<CompactPatchDescriptor> descriptors;
	vector<atom> const * atoms;
	vector<atom> outer_atoms;

	kdtree_atom atoms_tree;
	kdtree_atom core_atoms_tree;
	kdtree_atom outer_atoms_tree;

	float max_atm_radius;
	float min_atm_radius;
	vector<point3D> sphere_points;
	int n_sphere_points;

	Molecule(float patchRadius, float minCenterDist,
			float probeRadius, float resolution, string const & inname,
			string const & outname, string const & inname_radii,
			int maxOrder, bool no_hydrogen, bool no_hetatm,  bool receptor = true);
	virtual ~Molecule();
	void outputMoleculeDescriptors(string const & filename);
	void outputMoleculePQR(string const & filename);
	void outputMoleculePQR(string const & filename, DockingPose const & p);

	void calculateDescriptors(string const & outname, float minCenterDist, float patchRadius, int maxOrder, array3D const & interface, vector<CompactPatchDescriptor> & out_descriptors) {
		auto t_start = chrono::high_resolution_clock::now();
		cout << "Extracting surface patches\n";
		surface->extractPatchCenters(minCenterDist);
		cout << "Calculating Surface Descriptors\n";
		surface->calculateSurfaceDescriptors(patchRadius, maxOrder, interface, out_descriptors);
		auto t_ms = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t_start).count();
		cout << "Patch extraction and descriptor \ncalculation time:\t" << t_ms / 1000.0 << " seconds.\n"<<endl;
	}

//	void calculateInterfaceDescriptors(string const & outname, float minCenterDist, float patchRadius, int maxOrder, array3D const & interface, double threshold, vector<CompactPatchDescriptor> & out_descriptors) {
//		auto t_start = chrono::high_resolution_clock::now();
//		cout << "Extracting interface patches\n";
//		surface->extractInterfacePatchCenters(minCenterDist, interface);
//		cout << "Calculating interface descriptors\n";
//		surface->calculateSurfaceDescriptors(patchRadius, maxOrder, interface, out_descriptors);
//		auto t_ms = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t_start).count();
//		out_descriptors.erase(std::remove_if(out_descriptors.begin(), out_descriptors.end(),
//		                         [](CompactPatchDescriptor const & cpd)
//		                         { return  (not cpd.isInterface); }),
//				out_descriptors.end());
//		cout << "Patch extraction and descriptor \ncalculation time:\t" << t_ms / 1000.0 << " seconds.\n"<<endl;
//	}

	/**
	 * Returns list of coordinates on a sphere using the Golden-Section Spiral algorithm.
	 * @param n			number of points on the sphere
	 * @param points	output vector of generated points
	 */
	inline void generateSpherePoints(size_t n) {
		sphere_points.resize(n);
		double y, r, phi;
		double inc = M_PI * (3 - sqrt(5.0));
		double offset = 2.0 / n;
		for (size_t i = 0; i < n; ++i) {
		    y = i * offset - 1 + (offset / 2.0);
		    r = sqrt(1 - y * y);
		    phi = i * inc;
		    sphere_points[i] = point3D(cos(phi)*r, y, sin(phi)*r);
		}
	};
	/**
	 * Determines if the given point is buried inside the volume of the molecule.
	 * @param point		The input point
	 * @return		true if the input point is buried, false otherwise, i.e. if it is
	 * 				solvent accessible
	 */
	inline bool is_buried(point3D const & point) {
		vector<pair<size_t, float> > ret_matches;
		float searchRad = (this->surface->probeRadius + Molecule::max_atm_radius);
		float s_searchRad = searchRad * searchRad;
		size_t k = atoms_tree.radiusSearch(point, s_searchRad, ret_matches);
		for (size_t ii = 0; ii < k; ++ii) {
			atom const * nb = &this->pqrModel->atomsInModel[ret_matches[ii].first];
			float atm_rad = nb->radius + this->surface->probeRadius;
			if (floatCompare(nb->distance(point), atm_rad))
				continue;
			if (nb->distance(point) < atm_rad)
				return true;
		}
		return false;
	};
	inline double calculate_AccessibleSurfaceArea(vector<atom> const & atoms,
			vector<vector<point3D>> & per_atom_sa_points,
			vector<double> & per_atom_asa, size_t n_sphere_points = 96) {
		size_t num_atoms = atoms.size();
		if (num_atoms == 0)
			return -1;
		per_atom_sa_points.resize(num_atoms);
		per_atom_asa.resize(num_atoms);
#pragma omp critical
	{
		if (sphere_points.empty())
			generateSpherePoints(n_sphere_points);
	}
		double total_ASA = 0;

		double c = 4.0 * M_PI / n_sphere_points;
			for (size_t i = 0; i < num_atoms; ++i) {
				size_t n_accessible_pts = 0;
				double radius = atoms[i].radius + this->surface->probeRadius;
				point3D atomCenter(atoms[i].x, atoms[i].y, atoms[i].z);
				for (size_t j = 0; j < n_sphere_points; ++j) {
					point3D currentPoint = radius * sphere_points[j] + atomCenter;
					if (!is_buried(currentPoint)) {
						++n_accessible_pts;
						per_atom_sa_points[i].push_back(currentPoint);
					}
				}
				per_atom_asa[i] = c * n_accessible_pts * radius * radius;
				total_ASA += per_atom_asa[i];
			}

		return total_ASA;
	};
};


#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_MOLECULE_MOLECULE_H_ */
