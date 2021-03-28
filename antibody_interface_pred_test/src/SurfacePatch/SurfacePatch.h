/*
 * SurfacePatch.h
 *
 *  Created on: Jul 16, 2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_SURFACEPATCH_SURFACEPATCH_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_SURFACEPATCH_SURFACEPATCH_H_

#include "../Geometry/eig3/eig3.h"
#include "../Geometry/voxelGrid.h"
#include "../Geometry/voxel_offset.h"
#include "../Zernike/ZernikeDescriptor.h"
#include <list>

typedef struct SurfacePatch {
public:
	bool interface;
	size_t ID;	/**< The current Patch's ID */
	point3D COG;
	uint16_t patchRadius;
	voxelGrid const * molecularSurface;
	voxelGrid const * volumetricModel;
	array3D surface; /**< Pointer to the cubic grid containing the current patch.*/
	array3D potentials_pos;
	array3D potentials_neg;
	array3D hydrophobicity_pos, hydrophobicity_neg;

	array3D BLAM930101_pos;
	array3D BLAM930101_neg;
	array3D BIOV880101_pos;
	array3D BIOV880101_neg;
	array3D MAXF760101_;
	array3D TSAJ990101_;
	array3D NAKH920108_;
	array3D CEDJ970104_;
	array3D LIFS790101_;
	array3D MIYS990104_pos;
	array3D MIYS990104_neg;

	ZernikeDescriptor* surfaceD; 		/**< Zernike Descriptor of the surface of the current patch. */
	ZernikeDescriptor* hydrophobicity_posD; 		/**< Zernike Descriptor of the surface of the current patch. */
	ZernikeDescriptor* hydrophobicity_negD; 		/**< Zernike Descriptor of the surface of the current patch. */

	ZernikeDescriptor* potentials_posD; /**< Zernike Descriptor of the positive potentials distribution of the current patch. */
	ZernikeDescriptor* potentials_negD;	/**< Zernike Descriptor of the negative potentials distribution of the current patch. */

	ZernikeDescriptor* BLAM930101_posD;
	ZernikeDescriptor* BLAM930101_negD;
	ZernikeDescriptor* BIOV880101_posD;
	ZernikeDescriptor* BIOV880101_negD;
	ZernikeDescriptor* MAXF760101_D;
	ZernikeDescriptor* TSAJ990101_D;
	ZernikeDescriptor* NAKH920108_D;
	ZernikeDescriptor* CEDJ970104_D;
	ZernikeDescriptor* LIFS790101_D;
	ZernikeDescriptor* MIYS990104_posD;
	ZernikeDescriptor* MIYS990104_negD;

	point3D patchNormal;
	float curvature;
	float flatness;
	float resolution;
	voxel patchCenter;
	point3D translation;

	/**
	 * Constructor.
	 * @param id			Id of the current patch
	 * @param type 			patch type, i.e. concave, convex or flat
	 * @param patchRadius	radius of the patch
	 * @param center		center of the current patch
	 * @param resolution	resolution of the molecular surface
	 * @param molecularSurface		pointer to the original surface
	 * @param positivePotentials	pointer to the positive potentials grid of the current molecule
	 * @param negativePotentials	pointer to the negative potentials grid of the current molecule
	 * @param ptran		Translation to the original PDB coordinate syste
	 * @param maxOrder	maximal order for the Zernike Descriptor calculation
	 */

	SurfacePatch(size_t id, float patchRadius,
			voxel const & center, float resolution,
			voxelGrid const * molecularSurface,
			voxelGrid const * volumetricModel,
			array3D const & potentials,
			array3D const & hydrophobicity,
			point3D const & ptran,
			array3D const & interface,
			array3D const & _BLAM930101,
			array3D const & _BIOV880101,
			array3D const & _MAXF760101,
			array3D const & _TSAJ990101,
			array3D const & _NAKH920108,
			array3D const & _CEDJ970104,
			array3D const & _LIFS790101,
			array3D const & _MIYS990104);

	virtual ~SurfacePatch();
	/**
	 * Calculates the Zernike Descriptors of the surface of the current patch.
	 * If hydrophobicity and/or electrostatic potentials information is available,
	 * this method calculates the Zernike Descriptors for these distributions too.
	 */
	void calculateZernikeDescriptors(int order);
	/**
	 * Calculates the normal vector and curvature of the current patch.
	 * The two optional parameters are needed for the determination of the orientation of
	 * the estimated normal vector. The estimated normal is first scaled by the 'scale'
	 * parameter, then it is added to the current patch center. If the resulting voxel is
	 * outside the volumetric model of the molecule, then the normal orientation is correct.
	 * Otherwise, the normal vector is multiplied by -1. The normal has to be scaled, for us
	 * to be sure that we are reaching a voxel other than the current patch center, (the estimated
	 * normal has unitary norm).
	 *
	 */
	void normalEstimator();
    /**
     * Extraction operator
     * The overload of operator<< that takes a std::ostream& as the left hand argument.
     * Since this operator takes the user-defined type as the right argument (b in a@b),
	 * it must be implemented as non-members.
     */
	friend inline std::ostream & operator<<(std::ostream & os, SurfacePatch const & sp) {
		os.precision(std::numeric_limits<double>::digits10 + 2);
		os.setf(ios::scientific, ios::floatfield);
		os << "Patch ID: " << sp.ID
				<< ", radius: " << sp.patchRadius
				<< ", translation: " << sp.translation
				<< ", resolution: " << sp.resolution
				<< ", center: "	<< sp.patchCenter
				<< ", normal: " << sp.patchNormal
//				<< ", curvature: " << sp.curvature
				<< "\n";
		if (sp.surfaceD != NULL)
			os << "Surface Zernike Invariant:\n" << *sp.surfaceD << "\n";
		if (sp.potentials_posD != NULL)
			os << "Positive potentials Zernike Invariant:\n" << *sp.potentials_posD << "\n";
		if (sp.potentials_negD != NULL)
			os << "Negative potentials Zernike Invariant:\n" << *sp.potentials_negD << "\n";
		return os;
	};
	//------DEBUG--------
	/**
	 * Exports the patch surface to the PCD file format.
	 * @param filename		output file name
	 */
	void exportPatchPCDModel(string const & filename);
	void exportPatchPotentialsPCDModel(string const & filename);
	void exportPatchPPotentialsPCDModel(string const & filename);
	void exportPatchNPotentialsPCDModel(string const & filename);

	/**
	 * Exports the invariant(s) of the current surface patch
	 * to file.
	 * @param filename 		output file name
	 */
	void exportZernikeInvariants(string const & filename);
	/**
	 * Reconstructs the current patch surface from the precalculated
	 * Zernike Moments.
	 * @param threshold		threshold for the reconstructed grid value to
	 * 						be considered as 1.
	 */
	void reconstructPatch();
	/**
	 * Exports the reconstructed patch surface to the PCD file format.
	 * @param filename		output file name
	 */
	void exportReconstructedPatchPCDModel(string const & filename, double threshold);
private:
	int dim; 		/**< The dimension of the cubic grid containing the current patch.
			   	   	  *  The grid size is $dim^3$. */
	static list<voxel_offset> patchSphereOffsets; /**< list of offsets for the extraction of the patch */

	//-------debug-------
	array3D reconstructedPatch; /**< the cubic grid containing the current patch.*/
	/**
	 * Calulates the centroid (center of gravity) of the current patch.
	 * @return 	the centroid of the current patch
	 */
	inline void calculateCentroid() {
		COG = point3D(0, 0, 0);
		size_t numberOfVoxels = 0;
		voxel center(dim/2, dim/2, dim/2);
		for (auto const & offset : patchSphereOffsets) {
			voxel c = center + offset;
			if (fabs(surface[c.ix][c.iy][c.iz]) > 0) {
				COG.x += c.ix;
				COG.y += c.iy;
				COG.z += c.iz;
				++numberOfVoxels;
			}
		}
		if (numberOfVoxels != 0) /* avoid dividing by zero */
			COG /= numberOfVoxels;
	}; // calculateCentroid
	/**
	 * Calculates the correlation matrix for the current patch.
	 * This step is needed during the normal vector and curvature
	 * estimation process.
	 * @param correlation	the output correlation 3x3 matrix
	 */
	inline void calculateCorrelationMatrix(double correlation[3][3]) {
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				correlation[i][j] = 0;
		int numberOfVoxels = 0;
		voxel center(dim/2, dim/2, dim/2);
		calculateCentroid();
		for (auto const & offset : patchSphereOffsets) {
			voxel cv = center + offset;
			if (fabs(surface[cv.ix][cv.iy][cv.iz]) > 0) {
				point3D c(cv.ix, cv.iy, cv.iz);
				c -= COG;
				correlation[0][0] += c.x * c.x;
				correlation[0][1] += c.x * c.y;
				correlation[0][2] += c.x * c.z;
				correlation[1][0] += c.y * c.x;
				correlation[1][1] += c.y * c.y;
				correlation[1][2] += c.y * c.z;
				correlation[2][0] += c.z * c.x;
				correlation[2][1] += c.z * c.y;
				correlation[2][2] += c.z * c.z;
				++numberOfVoxels;
			}
		}
		if (numberOfVoxels != 0) { /* avoid dividing by zero */
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					correlation[i][j] /= numberOfVoxels;
		}
	} // calculateCorrelationMatrix

	/**
	 * This method clears all disconnected voxels that are not reachable from the
	 * center of the current patch.
	 */
	void clearDisconnectedVoxels();
} SurfacePatch;


#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_SURFACEPATCH_SURFACEPATCH_H_ */
