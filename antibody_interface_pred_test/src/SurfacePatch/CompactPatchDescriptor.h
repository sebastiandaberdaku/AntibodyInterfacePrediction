/*
 * CompactPatchDescriptor.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef COMPACTPATCHDESCRIPTOR_H_
#define COMPACTPATCHDESCRIPTOR_H_

#include "../utils/vectorDistance.h"
#include "SurfacePatch.h"
#include <assert.h>
#include <vector>

using namespace std;

/**
 * Struct containing the descriptors and other critical information of a SurfacePatch.
 */
typedef struct CompactPatchDescriptor {
	bool isInterface;
    size_t ID;      /**< The current Patch's ID */
	point3D center;

//  float curvature;
//	point3D patchNormal;
//	array1D surfaceI; 		 /**< Zernike Invariants of the surface of the current patch. */
//	array1D potentials_posI; /**< Zernike Invariants of the positive potentials distribution of the current patch. */
//	array1D potentials_negI; /**< Zernike Invariants of the negative potentials distribution of the current patch. */
//	array1D hydrophobicity_posI; /**< Zernike Invariants of the hydrophobicity distribution of the current patch. */
//	array1D hydrophobicity_negI; /**< Zernike Invariants of the hydrophobicity distribution of the current patch. */

	array1D BLAM930101_posI;
	array1D BLAM930101_negI;
	array1D BIOV880101_posI;
	array1D BIOV880101_negI;
	array1D MAXF760101_I;
	array1D TSAJ990101_I;
	array1D NAKH920108_I;
	array1D CEDJ970104_I;
	array1D LIFS790101_I;
	array1D MIYS990104_posI;
	array1D MIYS990104_negI;

	/**
	 * Default constructor.
	 */
	CompactPatchDescriptor();

	/**
	 * Copy constructor
	 * @param cpd	CompactPatchDescriptor struct to be copied
	 */
	CompactPatchDescriptor(CompactPatchDescriptor const & cpd);
	/**
	 * Copy assignment operator
	 */
	CompactPatchDescriptor & operator=(CompactPatchDescriptor const & cpd);
	/**
	 * This constructor creates the CompactPatchDescriptor from an input SurfacePatch
	 * @param sp	input SurfacePatch from which to create the CompatcPatchDescriptor
	 */
	CompactPatchDescriptor(SurfacePatch const & sp);

	/**
	 * Calculates the Pearson's correlation coefficient between the electrostatic
	 * potentials descriptor of the current patch and the input one. There are two
	 * descriptors for the electrostatic potential of a surface patch, one for the
	 * positive potentials and one for the negative potentials. The correlation
	 * is calculated matching the positive potentials descriptor of the first patch
	 * with the negative potentials descriptor of the second patch, and the negative
	 * potentials descriptor of the first patch with the positive potentials descriptor
	 * of the second patch.
	 * @param cpd	the input CompactPatchDescriptor
	 * @return		a value between +1 and −1 inclusive, where 1 is total positive
	 * 				correlation, 0 is no correlation, and −1 is total negative correlation
	 */
//	inline double potentialsPN_correlation(CompactPatchDescriptor const & cpd) const {
//		assert(this->potentials_posI.size() > 0
//				&& this->potentials_posI.size() == cpd.potentials_negI.size());
//		return correlationPearson(this->potentials_posI, cpd.potentials_negI);
//	};
//	inline double potentialsNP_correlation(CompactPatchDescriptor const & cpd) const {
//		assert(this->potentials_negI.size() > 0
//				&& this->potentials_negI.size() == cpd.potentials_posI.size());
//		return correlationPearson(this->potentials_negI, cpd.potentials_posI);
//	};
	/**
	 * Calculates the Pearson's correlation coefficient between the geometric surface
	 * descriptor of the current patch and the input one.
	 * @param cpd	the input CompactPatchDescriptor
	 * @return		a value between +1 and −1 inclusive, where 1 is total positive
	 * 				correlation, 0 is no correlation, and −1 is total negative correlation
	 */
//	inline double surface_correlation(CompactPatchDescriptor const & cpd) const {
//		assert(this->surfaceI.size() > 0
//				&& this->surfaceI.size() == cpd.surfaceI.size());
//		return correlationPearson(this->surfaceI, cpd.surfaceI);
//	};
    /**
     * Extraction operator
     * The overload of operator<< that takes a std::ostream& as the left hand argument.
     * Since this operator takes the user-defined type as the right argument (b in a@b),
	 * it must be implemented as non-members.
     */
	friend std::ostream & operator<<(std::ostream& os, CompactPatchDescriptor const & cpd) {
		os.precision(std::numeric_limits<double>::digits10 + 2);
//		os.precision(5);
//		os << cpd.ID << "\t" << cpd.center.x << "\t" << cpd.center.y << "\t" << cpd.center.z << "\t"
//				<< cpd.curvature << "\t" <<cpd.isInterface << "\t" << cpd.surfaceI.size();
		if (cpd.isInterface)
			os << "+1";
		else
			os << "-1";
		size_t c = 0;
//		for (auto const & d : cpd.surfaceI)
//			os << "\t" << ++c << ":" << d;
//		for (auto const & d : cpd.potentials_posI)
//			os << "\t" << ++c << ":" << d;
//		for (auto const & d : cpd.potentials_negI)
//			os << "\t" << ++c << ":" << d;
//		for (auto const & d : cpd.hydrophobicity_posI)
//			os << "\t" << ++c << ":" << d;
//		for (auto const & d : cpd.hydrophobicity_negI)
//			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.BLAM930101_posI)
			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.BLAM930101_negI)
			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.BIOV880101_posI)
			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.BIOV880101_negI)
			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.MAXF760101_I)
			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.TSAJ990101_I)
			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.NAKH920108_I)
			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.CEDJ970104_I)
			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.LIFS790101_I)
			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.MIYS990104_posI)
			os << "\t" << ++c << ":" << d;
		for (auto const & d : cpd.MIYS990104_negI)
			os << "\t" << ++c << ":" << d;

		return os;
	};


	static void to_ostream(std::ostream& os, CompactPatchDescriptor const & cpd, size_t order) {
		os.precision(std::numeric_limits<double>::digits10 + 2);
//		os.precision(5);
//		os << cpd.ID << "\t" << cpd.center.x << "\t" << cpd.center.y << "\t" << cpd.center.z << "\t"
//				<< cpd.curvature << "\t" <<cpd.isInterface << "\t" << cpd.surfaceI.size();
		if (cpd.isInterface)
			os << "+1";
		else
			os << "-1";

		size_t max;
		if (order % 2 == 0)
			max = (order + 2) * (order + 2) / 4;
		else
			max = (order + 1) * (order + 3) / 4;

		size_t c = 0;
//		for (auto const & d : cpd.surfaceI)
//			os << "\t" << ++c << ":" << d;
//		for (auto const & d : cpd.potentials_posI)
//			os << "\t" << ++c << ":" << d;
//		for (auto const & d : cpd.potentials_negI)
//			os << "\t" << ++c << ":" << d;
//		for (auto const & d : cpd.hydrophobicity_posI)
//			os << "\t" << ++c << ":" << d;
//		for (auto const & d : cpd.hydrophobicity_negI)
//			os << "\t" << ++c << ":" << d;
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.BLAM930101_posI[i];
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.BLAM930101_negI[i];
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.BIOV880101_posI[i];
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.BIOV880101_negI[i];
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.MAXF760101_I[i];
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.TSAJ990101_I[i];
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.NAKH920108_I[i];
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.CEDJ970104_I[i];
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.LIFS790101_I[i];
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.MIYS990104_posI[i];
		for (size_t i = 0; i < max; ++i)
			os << "\t" << ++c << ":" << cpd.MIYS990104_negI[i];
	};

	/**
	 * Equal?
	 */
	inline bool operator==(CompactPatchDescriptor const & rhs) const {
		return (this->ID == rhs.ID) && (this->center == rhs.center);
	};


} CompactPatchDescriptor;

namespace std {
template<>
struct hash<CompactPatchDescriptor> {
	size_t operator()(CompactPatchDescriptor const & cpd) const {
		size_t result = hash<size_t>()(cpd.ID);
		return result;
	}
};

template<>
struct equal_to<CompactPatchDescriptor> {
	bool operator()(CompactPatchDescriptor const & p, CompactPatchDescriptor const & q) const {
		return (p.ID == q.ID) && (p.center == q.center);
	}
};
}

#endif /* COMPACTPATCHDESCRIPTOR_H_ */
