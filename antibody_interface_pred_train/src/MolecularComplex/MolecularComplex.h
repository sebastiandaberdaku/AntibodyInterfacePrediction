/*
 * MolecularComplex.h
 *
 *  Created on: 03/nov/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_MOLECULARCOMPLEX_MOLECULARCOMPLEX_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_MOLECULARCOMPLEX_MOLECULARCOMPLEX_H_

#include "../ACE/AtomicContactEnergies.h"
#include "../ASP/AtomicSolvationParameters.h"
#include "../DockingMethods/DockingPose.h"
#include "../Molecule/Molecule.h"
#include <map>
#include <string.h>
#include <unordered_set>


struct residue{
	residue(char id, int n) : chain_ID(id), residue_number(n) { };
	char chain_ID; /**< An optional string which provides the chain ID of the atom. */
	int residue_number; /**< An integer which provides the residue index. */
};

namespace std {
template<>
struct hash<residue> {
	size_t operator()(residue const & p) const {
		size_t result = hash<char>()(p.chain_ID);
		size_t h1 = hash<int>()(p.residue_number);
		boost::hash_combine(result, h1);
		return result;
	}
};
template<>
struct equal_to<residue> {
	bool operator()(residue const & p, residue const & q) const {
		return (p.chain_ID == q.chain_ID) && (p.residue_number == q.residue_number);
	}
};	}


using namespace std;

class MolecularComplex {
public:
	Molecule const * receptor;
	array3D receptor_interface;

	unordered_set<residue> interface_receptor_residues;


	Molecule const * ligand;
	array3D ligand_interface;

	DockingPose const * dockingPose;

	MolecularComplex(Molecule const & receptor, Molecule const & ligand, float interface_distance);
	MolecularComplex(Molecule const & receptor, Molecule const & ligand, DockingPose const & p);
	virtual ~MolecularComplex();

	void outputSurfacePCDModel(string const & filename);
	bool calculateScores(double & hydrophobicity, double & bindingEnergy,
			double & CoulombicPotential, double & ACE);
	bool compenetrationAllowed();
	void outputReceptorInterfacePCDModel(string const & filename);
	void outputLigandInterfacePCDModel(string const & filename);

	inline bool is_rInterface(point3D const & p) const {
		/* Translate and discretize the coordinates */
		uint16_t cx = static_cast<uint16_t>((p.x + receptor->translation.x) * receptor->surface->resolution + 0.5);
		uint16_t cy = static_cast<uint16_t>((p.y + receptor->translation.y) * receptor->surface->resolution + 0.5);
		uint16_t cz = static_cast<uint16_t>((p.z + receptor->translation.z) * receptor->surface->resolution + 0.5);
		voxel c(cx, cy, cz);
		voxel nb;
		for (int i = 0; i < 26; ++i) {
			nb = c + voxel_offset::neighbours[i];
			if (receptor_interface[nb.ix][nb.iy][nb.iz] > 0)
				return true;
		}
		return false;
	};
	inline bool is_lInterface(point3D const & p) const {
		/* Translate and discretize the coordinates */
		uint16_t cx = static_cast<uint16_t>((p.x + ligand->translation.x) * ligand->surface->resolution + 0.5);
		uint16_t cy = static_cast<uint16_t>((p.y + ligand->translation.y) * ligand->surface->resolution + 0.5);
		uint16_t cz = static_cast<uint16_t>((p.z + ligand->translation.z) * ligand->surface->resolution + 0.5);
		voxel c(cx, cy, cz);
		voxel nb;
		for (int i = 0; i < 26; ++i) {
			nb = c + voxel_offset::neighbours[i];
			if (ligand_interface[nb.ix][nb.iy][nb.iz] > 0)
				return true;
		}
		return false;
	};
private:
	/**
	 * Applies a rotation and a translation to the all the atoms in the input vector.
	 * @param p		the DockingPose to be applied (rotation + translation)
	 */
	static inline void transformAtoms(vector<atom> & atoms, DockingPose const & p) {
		if (atoms.empty())
			return;
		for(auto & atm : atoms) {
			atm.transform(p.rotation, p.translation);
		}
	};
	/**
	 * Applies a rotation and a translation to the all the 3D points in the input vector.
	 * @param p		the DockingPose to be applied (rotation + translation)
	 */
	static inline void transformPoints(vector<point3D> & points, DockingPose const & p) {
		if (points.empty())
			return;
		vector<point3D>::iterator it;
		point3D temp;
		for (it = points.begin(); it !=  points.end(); ++it) {
			temp = point3D(it->x, it->y, it->z);
			temp = temp.transform(p.rotation);
			temp += p.translation;
			it->x = temp.x;
			it->y = temp.y;
			it->z = temp.z;
		}
	};
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
	 * @param ligand_atoms vector containing the ligand atoms transformed according
	 * 					the input docking pose
	 * @return true if the compenetration is acceptable, false otherwise
	 */
//	inline bool compenetrationAllowed(vector<atom> const & ligand_atoms) {
//		vector<atom>::const_iterator lig_atm;
//		for (lig_atm = ligand_atoms.begin(); lig_atm != ligand_atoms.end(); ++lig_atm) {
//			if (lig_atm->hasZeroRadius())
//				continue;
//			float searchRad = 2	* (receptor->surface->probeRadius + Molecule::max_atm_radius);
//			float ssRad = searchRad * searchRad;
//			vector<pair<size_t, float> > ret_matches;
//			size_t k = receptor->core_atoms_tree.radiusSearch((*lig_atm), ssRad, ret_matches);
//			for (size_t ii = 0; ii < k; ++ii) {
//				atom const * const rec_atm = &receptor->outer_atoms[ret_matches[ii].first];
//				if (rec_atm->hasZeroRadius())
//					continue;
//				float cutoff = (*lig_atm).radius + (*rec_atm).radius;
//				float s_cutoff = cutoff * cutoff;
//				float s_dist = (*lig_atm).s_distance((*rec_atm));
//				if (s_dist < s_cutoff)
//					return false;
//			}
//		}
//		return true;
//	};
//	/**
//	 * Calculate the hydrophobic energy term.
//	 * The receptor's atoms are inserted into a kdtree data structure;
//	 * For every ligand atom, we count the receptor atoms within a certain
//	 * distance, and, based on their charge we update three counters.
//	 * @param ligand_atoms vector containing the ligand atoms transformed according
//	 * 					the input docking pose
//	 * @return the hydrophobic energy term
//	 */
//	inline double getHydrophobicEnergy(vector<atom> const & ligand_atoms) {
//		int hh = 0; /**< number of contacts between non-polar (hydrophobic) atoms */
//		int pp = 0; /**< number of contacts between polar (hydrophilic) atoms */
//		int hp = 0; /**< number of contacts between non-polar (hydrophobic) and polar (hydrophilic) atoms */
//		float searchRad = 2	* (receptor->surface->probeRadius + Molecule::max_atm_radius);
//		float ssRad = searchRad * searchRad;
//		vector<atom>::const_iterator lig_atm;
//		for (lig_atm = ligand_atoms.begin(); lig_atm != ligand_atoms.end(); ++lig_atm) {
//			vector<pair<size_t, float> > ret_matches;
//			size_t k = receptor->outer_atoms_tree.radiusSearch((*lig_atm), ssRad, ret_matches);
//			for (size_t ii = 0; ii < k; ++ii) {
//				atom const * const rec_atm = &receptor->outer_atoms[ret_matches[ii].first];
//				float cutoff = 2.0 * receptor->surface->probeRadius + (*lig_atm).radius + (*rec_atm).radius;
//				float s_cutoff = cutoff * cutoff;
//				float s_dist = (*lig_atm).s_distance((*rec_atm));
//				if (s_dist >= s_cutoff)
//					continue;
//				if (floatCompare(lig_atm->charge, 0)) {
//					if (floatCompare(rec_atm->charge, 0))
//						++hh;
//					else
//						++hp;
//				} else {
//					if (floatCompare(rec_atm->charge, 0))
//						++hp;
//					else
//						++pp;
//				}
//			}
//		}
//		if (hh == 0)
//			return 0;
//		else
//			return hh / (double) (hh + pp + hp);
//	};
//	/**
//	 * Calculate the Coulombic potential energy term.
//	 * The receptor's atoms are inserted into a kdtree data structure;
//	 * For every ligand atom, we evaluate the electrostatic interaction
//	 * with the receptor atoms within a certain distance, and, based on
//	 * their charge we update the overall Coulombic Potential.
//	 * @param ligand_atoms vector containing the ligand atoms transformed according
//	 * 					the input docking pose
//	 * @return the Coulombic potential energy term
//	 */
//	inline double getCoulombicEnergy(vector<atom> const & ligand_atoms) {
//		double CoulombicPotential = 0;
//		float searchRad = 2	* (receptor->surface->probeRadius + Molecule::max_atm_radius);
//		float ssRad = searchRad * searchRad;
//		vector<atom>::const_iterator lig_atm;
//		for (lig_atm = ligand_atoms.begin(); lig_atm != ligand_atoms.end(); ++lig_atm) {
//			vector<pair<size_t, float> > ret_matches;
//			size_t k = receptor->outer_atoms_tree.radiusSearch((*lig_atm), ssRad, ret_matches);
//			for (size_t ii = 0; ii < k; ++ii) {
//				atom const * const rec_atm = &receptor->outer_atoms[ret_matches[ii].first];
//
//				float cutoff = 2.0 * receptor->surface->probeRadius + (*lig_atm).radius + (*rec_atm).radius;
//				float s_cutoff = cutoff * cutoff;
//				float s_dist = (*lig_atm).s_distance((*rec_atm));
//
//				if (s_dist >= s_cutoff)
//					continue;
//
//				double qiqj = lig_atm->charge * rec_atm->charge;
//				double d = sqrt(s_dist) + 1.5;
//				CoulombicPotential += qiqj / (d * d);
//			}
//		}
//		return CoulombicPotential;
//	};
//
//	/**
//	 * This method calculates the Atomic Contact Energy as decribed in:
//	 *
//	 * Determination of atomic desolvation energies from the structures of crystallized proteins.
//	 * Zhang C, Vasmatzis G, Cornette JL, DeLisi C
//	 * Department of Biomedical Engineering, Boston University, MA 02215, USA.
//	 * Journal of Molecular Biology [1997, 267(3):707-726]
//	 *
//	 * Instead of the rigid 6Å cutoff radius we dynamically determine if two atoms are in
//	 * contact. We consider only those atom-atom contacts that exclude an intervening water
//	 * molecule of radius 'probeRadius'.
//	 *
//	 * @param ligand_atoms vector containing the atoms of the ligand molecule
//	 * @return the value of the Atomic Contact Energy for the current pose
//	 */
//	inline float getAtomicContactEnergy(vector<atom> const & ligand_atoms) {
//		float ACE = 0;
//		float searchRad = 2	* (receptor->surface->probeRadius + Molecule::max_atm_radius);
//		float ssRad = searchRad * searchRad;
//		vector<atom>::const_iterator lig_atm;
//		for (lig_atm = ligand_atoms.begin(); lig_atm != ligand_atoms.end(); ++lig_atm) {
//			vector<pair<size_t, float> > ret_matches;
//			size_t k = receptor->outer_atoms_tree.radiusSearch((*lig_atm), ssRad, ret_matches);
//			for (size_t ii = 0; ii < k; ++ii) {
//				atom const * const rec_atm = &receptor->pqrModel->atomsInModel[ret_matches[ii].first];
//
//				float cutoff = 2.0 * receptor->surface->probeRadius + (*lig_atm).radius + (*rec_atm).radius;
//				float s_cutoff = cutoff * cutoff;
//				float s_dist = (*lig_atm).s_distance((*rec_atm));
//
//				if (s_dist >= s_cutoff)
//					continue;
//				atomType_t rec_type = AtomicContactEnergies::resolveAtomType(*rec_atm);
//				atomType_t lig_type = AtomicContactEnergies::resolveAtomType(*lig_atm);
//				float eij = AtomicContactEnergies::atom_contact_energies.find(make_pair(rec_type, lig_type))->second;
//				ACE += eij;
//			}
//		}
//		return ACE;
//	};
};


#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_MOLECULARCOMPLEX_MOLECULARCOMPLEX_H_ */
