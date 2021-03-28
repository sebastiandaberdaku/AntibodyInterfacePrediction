/*
 * atom.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef ATOM_H_
#define ATOM_H_

#include "../Geometry/point3D.h"
#include "../utils/numerical_limits.h"
#include "../utils/trim.h"
#include <iomanip>
#include <limits>
#include <map>
#include <math.h>
#include <ostream>

using namespace std;

/**
 * The ATOM records present the atomic coordinates for standard amino acids and
 * nucleotides. They also present the occupancy and temperature factor for each atom.
 * Non-polymer chemical coordinates use the HETATM record type. The element symbol is
 * always present on each ATOM record; charge is optional.
 */
#if defined TIGHT_PACKING
#pragma pack(push)  /* push current alignment to stack */
#pragma pack(1)     /* set alignment to 1 byte boundary */
#endif
typedef struct atom {
public:
	atom();
	atom(string const & field_name, int atom_number, string const & atom_name,
			string const & residue_name, char chain_ID, int residue_number,
			float x, float y, float z, float charge, float radius);
	atom(atom const & atm);
	atom(point3D const & pt);
	atom & operator=(atom const & atm);
	friend bool operator==(atom const & lhs, atom const & rhs);

	inline void transform(double const R[3][3], point3D const & T) {
		this->x += T.x;
		this->y += T.y;
		this->z += T.z;

		double x = R[0][0] * this->x + R[0][1] * this->y + R[0][2] * this->z;
		double y = R[1][0] * this->x + R[1][1] * this->y + R[1][2] * this->z;
		double z = R[2][0] * this->x + R[2][1] * this->y + R[2][2] * this->z;

		this->x = x;
		this->y = y;
		this->z = z;
	};

	bool is_hetero; /**< Specifies the type of PQR entry and should either be false for ATOM or true for HETATM. */
	int atom_number; /**< An integer which provides the atom index. */
	string atom_name; /**< A string which provides the atom name. */
	string residue_name; /**< A string which provides the residue name. */
	char chain_ID; /**< An optional string which provides the chain ID of the atom. */
	int residue_number; /**< An integer which provides the residue index. */
	float x, y, z; /**< 3 floats which provide the atomic coordiantes. */
	float charge; /**< A float which provides the atomic charge (in electrons). */
	float & occupancy = charge;
	float radius; /**< A float which provides the atomic radius (in Ångstrom). */
	float & tempFactor = radius;

	inline float s_distance(atom const & a) const {
		float dx = this->x - a.x;
		float dy = this->y - a.y;
		float dz = this->z - a.z;
		return (dx * dx + dy * dy + dz * dz);
	}
	inline float distance(atom const & a) const {
		return sqrt(this->s_distance(a));
	};
	inline float s_distance(point3D const & p) const {
		float dx = this->x - p.x;
		float dy = this->y - p.y;
		float dz = this->z - p.z;
		return (dx * dx + dy * dy + dz * dz);
	};
	inline float distance(point3D const & p) const {
		return sqrt(this->s_distance(p));
	};
	inline friend std::ostream & operator<<(std::ostream & os, atom const & atm) {
		if (atm.is_hetero)
			os << "HETATM";
		else
			os << "ATOM  ";
		os << setw(5) << atm.atom_number
				<< setw(5) << atm.atom_name
				<< setw(5) << atm.residue_name
				<< setw(2) << atm.chain_ID
				<< setw(4) << atm.residue_number
				<< setw(12) << std::fixed << setprecision(3) << atm.x
				<< setw(9) << atm.y
				<< setw(9) << atm.z
				<< setw(7) << setprecision(4) << atm.charge
				<< setw(7) << atm.radius;
		return os;
	};
	/**
	 * Adding atoms with zero radius and non-zero charge is a function of the PARSE forcefield,
	 * not the PDB2PQR server itself - if you try using AMBER or CHARMM you'll find that
	 * neither one gives hydrogens a zero radius.
	 *
	 * If you're interested in why PARSE assigns zero radii, you might want to check out the
	 * original paper: Sitkoff et al. 1994 J Phys Chem 98, 1978-88.
	 */
	inline bool hasZeroRadius() const {
		return floatCompare(this->radius, 0);
	};
	inline bool is_HETATM() const {
		return is_hetero;
	};
	inline bool is_Hydrogen() const {
		return (atom_name[0] == 'H');
	};
} atom;

#if defined TIGHT_PACKING
#pragma pack(pop)   /* restore original alignment from stack */
#endif

#endif /* ATOM_H_ */
