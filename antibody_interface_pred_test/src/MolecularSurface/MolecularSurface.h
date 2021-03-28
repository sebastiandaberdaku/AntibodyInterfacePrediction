/*
 * MolecularSurface.h
 *
 *  Created on: Nov 11, 2013
 *      Author: daberdaku
 */
/**
 * This header defines the MolecularSurface object, representing
 * the surface of a molecule, with the related methods for its
 * calculation.
 */
#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_MOLECULARSURFACE_MOLECULARSURFACE_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_MOLECULARSURFACE_MOLECULARSURFACE_H_

#include "../Geometry/point3D.h"
#include "../Geometry/PotentialGridDX.h"
#include "../Geometry/voxel_offset.h"
#include "../Geometry/voxelGrid.h"
#include "../Atom/atom.h"
#include "../SurfacePatch/CompactPatchDescriptor.h"
#include "../SurfacePatch/SurfacePatch.h"
#include "../hydrophobicity/hydrophobicity.h"
#include <list>
#include <queue>
#include <vector>
#include "../utils/hash.h"

/**
 * Dielectric constant (default 4.0)
 */
#define EPSILON 4.0

class MolecularSurface {
	friend class Molecule;
	friend class MolecularComplex;
public:
	voxelGrid *cpkModel, *surface;
	float resolution;
	array3D potentialsDX;
	std::vector<voxel> patchCenters;

	MolecularSurface(std::vector<atom> const & atoms, float probeRadius,
			float resolution, uint16_t length, uint16_t width, uint16_t height,
			point3D translation);
	virtual ~MolecularSurface();

	static void boundingBox(std::vector<atom> const & atomsInModel,
			float probeRadius, float patchRadius, float resolution, uint16_t & length,
			uint16_t & width, uint16_t & height, point3D & translation,
			float & max_atm_radius, float & min_atm_radius);
	void createCPKModel();
	void createCPKModel_old();

	void fillInternalCavities();
	void fillInternalCavities_old();

	void buildSurface();
	void buildSurface_old();

	void fastRegionGrowingEDT();
	void outputSurfacePCDModel(std::string const & filename);

	void calculateHydrophobicity(std::map<std::string, float> const & hydrophobicity);
	void calculateHydrophobicity(std::map<std::string, float> const & hydrophobicity, size_t window_size);

	void calculateAtomHydrophobicity(std::unordered_map<pair<string, string>, int> const & hydrophobicity);

	void outputHydrophobicity(std::string const & filename);

	void calculateCoulombicPotentials();
	void outputCoulombicPotentials(std::string const & filename);

	void calculateDXPotentials(std::string const & input_DX);
	void outputDXPotentials(std::string const & outname);

	void extractPatchCenters(float minDistance);
	void extractPatchCenters(float minDistance, array3D const & interface, vector<voxel> & patchCenters);

	void outputPatchCentersNormalsPCDModel(std::string const & filename, float patchRadius, float normalRadius);
	void outputPatchCentersPCDModel(std::string const & filename);

	void calculateSurfaceDescriptors(float patchRadius, int maxOrder, array3D const & interface, vector<CompactPatchDescriptor> & descriptors);

	void calculateSurfaceDescriptors(float patchRadius, int maxOrder, vector<voxel> const & patchCenters, array3D const & interface, vector<CompactPatchDescriptor> & descriptors);


	void extractInterfacePatchCenters(float minDistance, array3D const & interface);

	void extractInterfacePatchCenters(float minDistance, array3D const & interface, vector<voxel> & patchCenters);


	void outputSurfacePatchPCDModel(size_t id, float patchRadius, string const & outname);

	/**
	 * Simple molecular volume calculation method. Counts the
	 * voxels occupied by the molecule in its volumetric model,
	 * then divides the result by $res^3$ to get an estimate
	 * of the volume in cubic Å.
	 */
	inline double calculateSolventExcludedVolume() {
		assert(cpkModel != NULL);
		size_t v = 0;
		for (int i = 0; i < plength; ++i)
			for (int j = 0; j < pwidth; ++j)
				for (int k = 0; k < pheight; ++k)
					if (cpkModel->getVoxel(i, j, k))
						++v;
		return v * pow((double) resolution, -3.0);
	};
	double calculateAccessibleSurfaceArea(std::vector<double> & per_atom_asa, std::vector<std::vector<point3D>> & per_atom_sap, size_t n_sphere_points = 256);
	double calculateSolvationEnergy(std::map<std::string, float> const & asp, std::vector<double> const & per_atom_asa);
	double calculateSolvationEnergy(std::map<std::string, float> const & asp);
	double calculateHydrophobicEnergy(std::map<std::string, float> const & hydrophobicity, std::vector<double> const & per_atom_asa);

	array3D BLAM930101_;
	array3D BIOV880101_;
	array3D MAXF760101_;
	array3D TSAJ990101_;
	array3D NAKH920108_;
	array3D CEDJ970104_;
	array3D LIFS790101_;
	array3D MIYS990104_;
	void calculateHQI8();


private:
	point3D ptran;
	std::vector<atom> const * atoms;
	float probeRadius;
	uint16_t pheight, pwidth, plength;

	array3D hydrophobicityMap;
	array3D coulombicPotentials;
	PotentialGridDX * DXPotentials;

	vector<point3D> sphere_points;

	void floodFill3D(voxelGrid & grid, voxel const & startingVoxel);

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
	 * Simple method that tests if the current surface voxel is Solvent Accessible.
	 * @param v		the current voxel
	 * @return		true if v is Solvent Accessible, false otherwise
	 */
	inline bool is_SA(voxel const & v) {
		voxel nb;
		for (int i = 0; i < 26; ++i) {
			nb = v + voxel_offset::neighbours[i];
			if (!cpkModel->getVoxel(nb))
				return true;
		}
		return false;
	};

	inline bool is_SA(point3D const & p) {
		/* Translate and discretize the coordinates */
		uint16_t cx = static_cast<uint16_t>((p.x + ptran.x) * resolution + 0.5);
		uint16_t cy = static_cast<uint16_t>((p.y + ptran.y) * resolution + 0.5);
		uint16_t cz = static_cast<uint16_t>((p.z + ptran.z) * resolution + 0.5);
		voxel c(cx, cy, cz);
		return is_SA(c);
	};
};
#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_MOLECULARSURFACE_MOLECULARSURFACE_H_ */
