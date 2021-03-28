/*
 * AtomicContactEnergies.h
 *
 *  Created on: 12/nov/2014
 *      Author: sebastian
 */
/**
 * Determination of atomic desolvation energies from the structures of crystallized proteins.
 * Zhang C, Vasmatzis G, Cornette JL, DeLisi C
 * Department of Biomedical Engineering, Boston University, MA 02215, USA.
 * Journal of Molecular Biology [1997, 267(3):707-726]
 */
#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_ACE_ATOMICCONTACTENERGIES_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_ACE_ATOMICCONTACTENERGIES_H_

#include "../utils/str2int.h"
#include "../utils/trim.h"
#include <map>

#include "../Atom/atom.h"
using namespace std;
/**
 * Only heavy atoms were considered in our computations and all hydrogen atoms were
 * combined with the heavy atoms to which they are bound. 18 heavy-atom types were
 * distinguished. The 'OTHER_ATM' atom type was introduced in this work to distinguish
 * possible unknown atom types.
 */
typedef enum atomType_t{
	N, CA, C, O, GCA, CB, KNZ, KCD, DOD, RNH, NND, RNE, SOG, HNE, YCZ, FCZ, LCD, CSG, OTHER_ATM
} atomType_t;

typedef struct AtomicContactEnergies {

	static atomType_t resolveAtomType(atom const & atm);
	static const map<pair<atomType_t, atomType_t>, float> atom_contact_energies;
	static const map<pair<atomType_t, atomType_t>, float> atom_contact_energies_other;

} AtomicContactEnergies;

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_ACE_ATOMICCONTACTENERGIES_H_ */
