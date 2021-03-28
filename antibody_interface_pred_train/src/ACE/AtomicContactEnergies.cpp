/*
 * AtomicContactEnergies.cpp
 *
 *  Created on: 01/mar/2015
 *      Author: sebastian
 */

#include "AtomicContactEnergies.h"

/**
 * Relative contact energy eij defined as follows:
 * 			eij = Eij + E00 - Ei0 - Ej0
 * where Eij is the absolute contact energy between atoms of
 * type i and j (0 - solvent atoms).
 *
 * By definition, eij is the energy change on forming an i-j pair
 * (solute-solute) pair and a 0-0 pair (solvent-solvent pair)
 * from an i-0 and a j-0 pair.
 *
 * In practice, eij is more appropriate than e'ij for calculating
 * the total contact energy of protein conformations because only
 * the number of atom-atom contacts, not including atom-solvent,
 * is used in computing that part of Ec attributable to protein
 * conformation.
 */
const map<pair<atomType_t, atomType_t>, float> AtomicContactEnergies::atom_contact_energies = {
		{ make_pair(N,	  N),	-0.724 },	{ make_pair(N,  CA),	-0.903 },	{ make_pair(N,	 C),	-0.722 },
		{ make_pair(N,	  O),	-0.322 },	{ make_pair(N, GCA),	-0.331 },	{ make_pair(N,	CB),	-0.603 },
		{ make_pair(N,	KNZ),	 1.147 },	{ make_pair(N, KCD),	 0.937 },	{ make_pair(N, DOD),	 0.752 },
		{ make_pair(N,	RNH),	 0.502 },	{ make_pair(N, NND),	 0.405 },	{ make_pair(N, RNE),	 0.495 },
		{ make_pair(N,	SOG),	-0.116 },	{ make_pair(N, HNE),	-0.092 },	{ make_pair(N, YCZ),	-0.412 },
		{ make_pair(N,	FCZ),	-0.499 },	{ make_pair(N, LCD),	-1.005 },	{ make_pair(N, CSG), 	-2.060 },

		{ make_pair(CA,	  N),	-0.903 },	{ make_pair(CA,	 CA),	-0.842 },	{ make_pair(CA,	  C),	-0.850 },
		{ make_pair(CA,	  O),	-0.332 },	{ make_pair(CA,	GCA),	-0.365 },	{ make_pair(CA,	 CB),	-0.531 },
		{ make_pair(CA,	KNZ),	 1.139 },	{ make_pair(CA, KCD),	 0.918 },	{ make_pair(CA, DOD),	 0.852 },
		{ make_pair(CA,	RNH),	 0.452 },	{ make_pair(CA, NND),	 0.445 },	{ make_pair(CA, RNE),	 0.436 },
		{ make_pair(CA,	SOG),	-0.044 },	{ make_pair(CA, HNE),	-0.165 },	{ make_pair(CA, YCZ),	-0.539 },
		{ make_pair(CA,	FCZ),	-0.590 },	{ make_pair(CA, LCD),	-1.123 },	{ make_pair(CA, CSG), 	-2.028 },

		{ make_pair(C,	  N),	-0.722 },	{ make_pair(C,  CA),	-0.850 },	{ make_pair(C,	 C),	-0.704 },
		{ make_pair(C,	  O),	-0.371 },	{ make_pair(C, GCA),	-0.308 },	{ make_pair(C,  CB),	-0.499 },
		{ make_pair(C,	KNZ),	 1.140 },	{ make_pair(C, KCD),	 0.846 },	{ make_pair(C, DOD),	 0.788 },
		{ make_pair(C,	RNH),	 0.490 },	{ make_pair(C, NND),	 0.408 },	{ make_pair(C, RNE),	 0.418 },
		{ make_pair(C,	SOG),	-0.045 },	{ make_pair(C, HNE),	-0.089 },	{ make_pair(C, YCZ),	-0.415 },
		{ make_pair(C,	FCZ),	-0.478 },	{ make_pair(C, LCD),	-0.963 },	{ make_pair(C, CSG), 	-2.033 },

		{ make_pair(O,	  N),	-0.322 },	{ make_pair(O,	CA),	-0.332 }, 	{ make_pair(O,	 C),	-0.371 },
		{ make_pair(O,	  O),	-0.016 },	{ make_pair(O, GCA),	 0.205 },	{ make_pair(O,  CB),	-0.059 },
		{ make_pair(O,	KNZ),	 1.414 },	{ make_pair(O, KCD),	 1.075 },	{ make_pair(O, DOD),	 1.152 },
		{ make_pair(O,	RNH),	 0.831 },	{ make_pair(O, NND),	 0.758 },	{ make_pair(O, RNE),	 0.649 },
		{ make_pair(O,	SOG),	 0.383 },	{ make_pair(O, HNE),	 0.241 },	{ make_pair(O, YCZ),	-0.073 },
		{ make_pair(O,	FCZ),	-0.166 },	{ make_pair(O, LCD),	-0.650 },	{ make_pair(O, CSG), 	-1.650 },

		{ make_pair(GCA,   N),	-0.331 },	{ make_pair(GCA,  CA),	-0.365 },	{ make_pair(GCA,   C),	-0.308 },
		{ make_pair(GCA,   O),	 0.205 },	{ make_pair(GCA, GCA),	 0.182 },	{ make_pair(GCA,  CB),	-0.009 },
		{ make_pair(GCA, KNZ),	 1.502 },	{ make_pair(GCA, KCD),	 1.385 },	{ make_pair(GCA, DOD),	 1.192 },
		{ make_pair(GCA, RNH),	 0.912 },	{ make_pair(GCA, NND),	 0.872 },	{ make_pair(GCA, RNE),	 0.943 },
		{ make_pair(GCA, SOG),	 0.434 },	{ make_pair(GCA, HNE),	 0.392 },	{ make_pair(GCA, YCZ),	-0.062 },
		{ make_pair(GCA, FCZ),	-0.021 },	{ make_pair(GCA, LCD),	-0.342 },	{ make_pair(GCA, CSG), 	-1.212 },

		{ make_pair(CB,	  N),	-0.603 },	{ make_pair(CB,	 CA),	-0.531 },	{ make_pair(CB,   C),	-0.499 },
		{ make_pair(CB,   O),	-0.059 },	{ make_pair(CB, GCA),	-0.009 },	{ make_pair(CB,  CB),	-0.469 },
		{ make_pair(CB, KNZ),	 1.310 },	{ make_pair(CB, KCD),	 1.010 },	{ make_pair(CB, DOD),	 0.859 },
		{ make_pair(CB, RNH),	 0.671 },	{ make_pair(CB, NND),	 0.591 },	{ make_pair(CB, RNE),	 0.565 },
		{ make_pair(CB, SOG),	 0.122 },	{ make_pair(CB, HNE),	 0.000 },	{ make_pair(CB, YCZ),	-0.438 },
		{ make_pair(CB, FCZ),	-0.533 },	{ make_pair(CB, LCD),	-1.034 },	{ make_pair(CB, CSG), 	-1.800 },

		{ make_pair(KNZ,   N),	 1.147 },	{ make_pair(KNZ,  CA),	 1.139 },	{ make_pair(KNZ,   C),	 1.140 },
		{ make_pair(KNZ,   O),	 1.414 },	{ make_pair(KNZ, GCA),	 1.502 },	{ make_pair(KNZ,  CB),	 1.310 },
		{ make_pair(KNZ, KNZ),	 3.018 },	{ make_pair(KNZ, KCD),	 2.911 },	{ make_pair(KNZ, DOD),	 1.157 },
		{ make_pair(KNZ, RNH),	 2.635 },	{ make_pair(KNZ, NND),	 1.848 },	{ make_pair(KNZ, RNE),	 2.699 },
		{ make_pair(KNZ, SOG),	 1.587 },	{ make_pair(KNZ, HNE),	 1.557 },	{ make_pair(KNZ, YCZ),	 1.004 },
		{ make_pair(KNZ, FCZ),	 1.340 },	{ make_pair(KNZ, LCD),	 1.252 },	{ make_pair(KNZ, CSG), 	 0.700 },

		{ make_pair(KCD,   N),	 0.937 },	{ make_pair(KCD,  CA),	 0.918 },	{ make_pair(KCD,   C),	 0.846 },
		{ make_pair(KCD,   O),	 1.075 },	{ make_pair(KCD, GCA),	 1.385 },	{ make_pair(KCD,  CB),	 1.010 },
		{ make_pair(KCD, KNZ),	 2.911 },	{ make_pair(KCD, KCD),	 2.811 },	{ make_pair(KCD, DOD),	 0.989 },
		{ make_pair(KCD, RNH),	 2.439 },	{ make_pair(KCD, NND),	 1.622 },	{ make_pair(KCD, RNE),	 2.525 },
		{ make_pair(KCD, SOG),	 1.361 },	{ make_pair(KCD, HNE),	 1.277 },	{ make_pair(KCD, YCZ),	 0.685 },
		{ make_pair(KCD, FCZ),	 0.938 },	{ make_pair(KCD, LCD),	 0.814 },	{ make_pair(KCD, CSG), 	 0.411 },

		{ make_pair(DOD,   N),	 0.752 },	{ make_pair(DOD,  CA),	 0.852 },	{ make_pair(DOD,   C),	 0.788 },
		{ make_pair(DOD,   O),	 1.152 },	{ make_pair(DOD, GCA),	 1.192 },	{ make_pair(DOD,  CB),	 0.859 },
		{ make_pair(DOD, KNZ),	 1.157 },	{ make_pair(DOD, KCD),	 0.989 },	{ make_pair(DOD, DOD),	-1.978 },
		{ make_pair(DOD, RNH),	 0.695 },	{ make_pair(DOD, NND),	 1.386 },	{ make_pair(DOD, RNE),	 0.641 },
		{ make_pair(DOD, SOG),	 1.002 },	{ make_pair(DOD, HNE),	 0.642 },	{ make_pair(DOD, YCZ),	 0.618 },
		{ make_pair(DOD, FCZ),	 0.894 },	{ make_pair(DOD, LCD),	 0.792 },	{ make_pair(DOD, CSG), 	-0.029 },

		{ make_pair(RNH,   N),	 0.502 },	{ make_pair(RNH,  CA),	 0.452 },	{ make_pair(RNH,   C),	 0.490 },
		{ make_pair(RNH,   O),	 0.831 },	{ make_pair(RNH, GCA),	 0.912 },	{ make_pair(RNH,  CB),	 0.671 },
		{ make_pair(RNH, KNZ),	 2.635 },	{ make_pair(RNH, KCD),	 2.439 },	{ make_pair(RNH, DOD),	 0.695 },
		{ make_pair(RNH, RNH),	 1.589 },	{ make_pair(RNH, NND),	 1.395 },	{ make_pair(RNH, RNE),	 1.515 },
		{ make_pair(RNH, SOG),	 1.022 },	{ make_pair(RNH, HNE),	 0.948 },	{ make_pair(RNH, YCZ),	 0.457 },
		{ make_pair(RNH, FCZ),	 0.707 },	{ make_pair(RNH, LCD),	 0.586 },	{ make_pair(RNH, CSG), 	 0.021 },

		{ make_pair(NND,   N),	 0.405 },	{ make_pair(NND,  CA),	 0.445 },	{ make_pair(NND,   C),	 0.408 },
		{ make_pair(NND,   O),	 0.758 },	{ make_pair(NND, GCA),	 0.872 },	{ make_pair(NND,  CB),	 0.591 },
		{ make_pair(NND, KNZ),	 1.848 },	{ make_pair(NND, KCD),	 1.622 },	{ make_pair(NND, DOD),	 1.386 },
		{ make_pair(NND, RNH),	 1.395 }, 	{ make_pair(NND, NND),	 1.300 },	{ make_pair(NND, RNE),	 1.283 },
		{ make_pair(NND, SOG),	 0.931 },	{ make_pair(NND, HNE),	 0.829 },	{ make_pair(NND, YCZ),	 0.484 },
		{ make_pair(NND, FCZ),	 0.633 },	{ make_pair(NND, LCD),	 0.389 },	{ make_pair(NND, CSG), 	-0.046 },

		{ make_pair(RNE,   N),	 0.495 },	{ make_pair(RNE,  CA),	 0.436 },	{ make_pair(RNE,   C),	 0.418 },
		{ make_pair(RNE,   O),	 0.649 },	{ make_pair(RNE, GCA),	 0.943 },	{ make_pair(RNE,  CB),	 0.565 },
		{ make_pair(RNE, KNZ),	 2.699 },	{ make_pair(RNE, KCD),	 2.525 },	{ make_pair(RNE, DOD),	 0.641 },
		{ make_pair(RNE, RNH),	 1.515 },	{ make_pair(RNE, NND),	 1.283 }, 	{ make_pair(RNE, RNE),	 1.309 },
		{ make_pair(RNE, SOG),	 0.921 },	{ make_pair(RNE, HNE),	 0.777 },	{ make_pair(RNE, YCZ),	 0.265 },
		{ make_pair(RNE, FCZ),	 0.484 },	{ make_pair(RNE, LCD),	 0.274 },	{ make_pair(RNE, CSG), 	-0.253 },

		{ make_pair(SOG,   N),	-0.116 },	{ make_pair(SOG,  CA),	-0.044 },	{ make_pair(SOG,   C),	-0.045 },
		{ make_pair(SOG,   O),	 0.383 },	{ make_pair(SOG, GCA),	 0.434 },	{ make_pair(SOG,  CB),	 0.122 },
		{ make_pair(SOG, KNZ),	 1.587 },	{ make_pair(SOG, KCD),	 1.361 },	{ make_pair(SOG, DOD),	 1.002 },
		{ make_pair(SOG, RNH),	 1.022 },	{ make_pair(SOG, NND),	 0.931 },	{ make_pair(SOG, RNE),	 0.921 },
		{ make_pair(SOG, SOG),	 0.559 },	{ make_pair(SOG, HNE),	 0.319 },	{ make_pair(SOG, YCZ),	 0.301 },
		{ make_pair(SOG, FCZ),	 0.212 },	{ make_pair(SOG, LCD),	-0.155 },	{ make_pair(SOG, CSG), 	-0.949 },

		{ make_pair(HNE,   N),	-0.092 },	{ make_pair(HNE,  CA),	-0.165 },	{ make_pair(HNE,   C),	-0.089 },
		{ make_pair(HNE,   O),	 0.241 },	{ make_pair(HNE, GCA),	 0.392 },	{ make_pair(HNE,  CB),	 0.000 },
		{ make_pair(HNE, KNZ),	 1.557 },	{ make_pair(HNE, KCD),	 1.277 },	{ make_pair(HNE, DOD),	 0.642 },
		{ make_pair(HNE, RNH),	 0.948 },	{ make_pair(HNE, NND),	 0.829 },	{ make_pair(HNE, RNE),	 0.777 },
		{ make_pair(HNE, SOG),	 0.319 },	{ make_pair(HNE, HNE),	-0.301 },	{ make_pair(HNE, YCZ),	-0.159 },
		{ make_pair(HNE, FCZ),	-0.132 },	{ make_pair(HNE, LCD),	-0.338 },	{ make_pair(HNE, CSG), 	-1.324 },

		{ make_pair(YCZ,   N),	-0.412 },	{ make_pair(YCZ,  CA),	-0.539 },	{ make_pair(YCZ,   C),	-0.415 },
		{ make_pair(YCZ,   O),	-0.073 },	{ make_pair(YCZ, GCA),	-0.062 },	{ make_pair(YCZ,  CB),	-0.438 },
		{ make_pair(YCZ, KNZ),	 1.004 },	{ make_pair(YCZ, KCD),	 0.685 },	{ make_pair(YCZ, DOD),	 0.618 },
		{ make_pair(YCZ, RNH),	 0.457 },	{ make_pair(YCZ, NND),	 0.484 },	{ make_pair(YCZ, RNE),	 0.265 },
		{ make_pair(YCZ, SOG),	 0.301 },	{ make_pair(YCZ, HNE),	-0.159 },	{ make_pair(YCZ, YCZ),	-0.314 },
		{ make_pair(YCZ, FCZ),	-0.478 },	{ make_pair(YCZ, LCD),	-0.964 },	{ make_pair(YCZ, CSG), 	-1.410 },

		{ make_pair(FCZ,   N),	-0.499 },	{ make_pair(FCZ,  CA),	-0.590 },	{ make_pair(FCZ,   C),	-0.478 },
		{ make_pair(FCZ,   O),	-0.166 },	{ make_pair(FCZ, GCA),	-0.021 },	{ make_pair(FCZ,  CB),	-0.533 },
		{ make_pair(FCZ, KNZ),	 1.340 },	{ make_pair(FCZ, KCD),	 0.938 },	{ make_pair(FCZ, DOD),	 0.894 },
		{ make_pair(FCZ, RNH),	 0.707 },	{ make_pair(FCZ, NND),	 0.633 },	{ make_pair(FCZ, RNE),	 0.484 },
		{ make_pair(FCZ, SOG),	 0.212 },	{ make_pair(FCZ, HNE),	-0.132 },	{ make_pair(FCZ, YCZ),	-0.478 },
		{ make_pair(FCZ, FCZ),	-0.687 },	{ make_pair(FCZ, LCD),	-1.240 },	{ make_pair(FCZ, CSG), 	-1.784 },

		{ make_pair(LCD,   N),	-1.005 },	{ make_pair(LCD,  CA),	-1.123 },	{ make_pair(LCD,   C),	-0.963 },
		{ make_pair(LCD,   O),	-0.650 },	{ make_pair(LCD, GCA),	-0.342 },	{ make_pair(LCD,  CB),	-1.034 },
		{ make_pair(LCD, KNZ),	 1.252 },	{ make_pair(LCD, KCD),	 0.814 },	{ make_pair(LCD, DOD),	 0.792 },
		{ make_pair(LCD, RNH),	 0.586 },	{ make_pair(LCD, NND),	 0.389 },	{ make_pair(LCD, RNE),	 0.274 },
		{ make_pair(LCD, SOG),	-0.155 },	{ make_pair(LCD, HNE),	-0.338 },	{ make_pair(LCD, YCZ),	-0.964 },
		{ make_pair(LCD, FCZ),	-1.240 },	{ make_pair(LCD, LCD),	-1.873 },	{ make_pair(LCD, CSG), 	-2.402 },

		{ make_pair(CSG,   N), 	-2.060 },	{ make_pair(CSG,  CA), 	-2.028 },	{ make_pair(CSG,   C), 	-2.033 },
		{ make_pair(CSG,   O), 	-1.650 },	{ make_pair(CSG, GCA), 	-1.212 },	{ make_pair(CSG,  CB), 	-1.800 },
		{ make_pair(CSG, KNZ), 	 0.700 },	{ make_pair(CSG, KCD), 	 0.411 },	{ make_pair(CSG, DOD), 	-0.029 },
		{ make_pair(CSG, RNH), 	 0.021 },	{ make_pair(CSG, NND), 	-0.046 },	{ make_pair(CSG, RNE), 	-0.253 },
		{ make_pair(CSG, SOG), 	-0.949 },	{ make_pair(CSG, HNE), 	-1.324 },	{ make_pair(CSG, YCZ), 	-1.410 },
		{ make_pair(CSG, FCZ), 	-1.784 },	{ make_pair(CSG, LCD), 	-2.402 },	{ make_pair(CSG, CSG), 	-3.742 }

};

/**
 * Relative contact energy e'ij defined as follows:
 * 			e'ij = Eij - (Eii + Ejj)/2
 * where Eij is the absolute contact energy between atoms of
 * type i and j (0 - solvent atoms).
 *
 * By definition, e'ij is half of the energy change due to forming
 * two i-j pairs from an i-i and a j-j pair.
 *
 * The expression for e'ij is more common than eij in lattice theories.
 */
const map<pair<atomType_t, atomType_t>, float> AtomicContactEnergies::atom_contact_energies_other = {
		{ make_pair(N,	 N),	 0.000 },	{ make_pair(N,  CA),	-0.119 },	{ make_pair(N,   C),	-0.008 },
		{ make_pair(N,	 O),	 0.048 },	{ make_pair(N, GCA),	-0.059 },	{ make_pair(N,  CB),	-0.007 },
		{ make_pair(N, KNZ),	 0.000 },	{ make_pair(N, KCD),	-0.106 },	{ make_pair(N, DOD),	 0.126 },
		{ make_pair(N, RNH),	 0.069 },	{ make_pair(N, NND),	 0.117 },	{ make_pair(N, RNE),	 0.203 },
		{ make_pair(N, SOG),	-0.033 },	{ make_pair(N, HNE),	 0.421 },	{ make_pair(N, YCZ),	 0.107 },
		{ make_pair(N, FCZ),	 0.207 },	{ make_pair(N, LCD),	 0.293 },	{ make_pair(N, CSG),	 0.173 },

		{ make_pair(CA,   N),	-0.119 },	{ make_pair(CA,  CA),	 0.000 },	{ make_pair(CA,   C),	-0.077 },
		{ make_pair(CA,   O),	 0.097 },	{ make_pair(CA, GCA),	-0.034 },	{ make_pair(CA,  CB),	 0.124 },
		{ make_pair(CA, KNZ),	 0.051 },	{ make_pair(CA, KCD),	-0.066 },	{ make_pair(CA, DOD),	 0.284 },
		{ make_pair(CA, RNH), 	 0.079 },	{ make_pair(CA, NND), 	 0.216 },	{ make_pair(CA, RNE),	 0.203 },
		{ make_pair(CA, SOG),	 0.097 },	{ make_pair(CA, HNE), 	 0.406 },	{ make_pair(CA, YCZ),	 0.039 },
		{ make_pair(CA, FCZ),	 0.175 },	{ make_pair(CA, LCD),	 0.235 },	{ make_pair(CA, CSG),	 0.264 },

		{ make_pair(C,   N),	-0.008 },	{ make_pair(C,  CA),	-0.077 },	{ make_pair(C,   C),	 0.000 },
		{ make_pair(C,   O),	-0.011 },	{ make_pair(C, GCA),	-0.047 },	{ make_pair(C,  CB),	 0.088 },
		{ make_pair(C, KNZ),	-0.017 },	{ make_pair(C, KCD), 	-0.207 },	{ make_pair(C, DOD), 	 0.151 },
		{ make_pair(C, RNH),     0.047 },	{ make_pair(C, NND),  	 0.110 },	{ make_pair(C, RNE), 	 0.116 },
		{ make_pair(C, SOG),	 0.028 },	{ make_pair(C, HNE),	 0.413 },	{ make_pair(C, YCZ), 	 0.094 },
		{ make_pair(C, FCZ),	 0.218 },	{ make_pair(C, LCD),  	 0.325 },	{ make_pair(C, CSG),	 0.190 },

		{ make_pair(O,   N), 	 0.048 },	{ make_pair(O,  CA),	 0.097 },	{ make_pair(O,   C),	-0.011 },
		{ make_pair(O,   O),	 0.000 },	{ make_pair(O, GCA),	 0.122 },	{ make_pair(O,  CB),	 0.184 },
		{ make_pair(O, KNZ),	-0.088 },	{ make_pair(O, KCD),	-0.323 },	{ make_pair(O, DOD),	 0.171 },
		{ make_pair(O, RNH),	 0.045 },	{ make_pair(O, NND),	 0.115 },	{ make_pair(O, RNE),	 0.002 },
		{ make_pair(O, SOG),	 0.111 },	{ make_pair(O, HNE),	 0.399 },	{ make_pair(O, YCZ),	 0.092 },
		{ make_pair(O, FCZ),	 0.185 },	{ make_pair(O, LCD),	 0.294 },	{ make_pair(O, CSG),	 0.228 },

		{ make_pair(GCA,   N),	-0.059 },	{ make_pair(GCA,  CA),	-0.034 },	{ make_pair(GCA,   C),	-0.047 },
		{ make_pair(GCA,   O),	 0.122 },	{ make_pair(GCA, GCA),	 0.000 },	{ make_pair(GCA,  CB),	 0.135 },
		{ make_pair(GCA, KNZ),	-0.098 },	{ make_pair(GCA, KCD),	-0.111 },	{ make_pair(GCA, DOD),	 0.112 },
		{ make_pair(GCA, RNH),	 0.027 },	{ make_pair(GCA, NND),	 0.131 },	{ make_pair(GCA, RNE),	 0.197 },
		{ make_pair(GCA, SOG),	 0.064 },	{ make_pair(GCA, HNE),	 0.451 },	{ make_pair(GCA, YCZ),	 0.004 },
		{ make_pair(GCA, FCZ),	 0.232 },	{ make_pair(GCA, LCD),	 0.504 },	{ make_pair(GCA, CSG),	 0.569 },

		{ make_pair(CB,   N),	-0.007 },	{ make_pair(CB,  CA),	 0.124 },	{ make_pair(CB,   C),	 0.088 },
		{ make_pair(CB,   O),	 0.184 },	{ make_pair(CB, GCA),	 0.135 },	{ make_pair(CB,  CB),	 0.000 },
		{ make_pair(CB, KNZ),	 0.036 },	{ make_pair(CB, KCD),	-0.161 },	{ make_pair(CB, DOD),	 0.105 },
		{ make_pair(CB, RNH),	 0.111 },	{ make_pair(CB, NND),	 0.175 },	{ make_pair(CB, RNE),	 0.145 },
		{ make_pair(CB, SOG),	 0.077 },	{ make_pair(CB, HNE),	 0.385 },	{ make_pair(CB, YCZ),	-0.046 },
		{ make_pair(CB, FCZ),	 0.045 },	{ make_pair(CB, LCD),	 0.137 },	{ make_pair(CB, CSG),	 0.306 },

		{ make_pair(KNZ,   N),	 0.000 },	{ make_pair(KNZ,  CA),	 0.051 },	{ make_pair(KNZ,   C),	-0.017 },
		{ make_pair(KNZ,   O),	-0.088 },	{ make_pair(KNZ, GCA),	-0.098 },	{ make_pair(KNZ,  CB),	 0.036 },
		{ make_pair(KNZ, KNZ),	 0.000 },	{ make_pair(KNZ, KCD),	-0.003 },	{ make_pair(KNZ, DOD),	-1.341 },
		{ make_pair(KNZ, RNH),	 0.331 },	{ make_pair(KNZ, NND),	-0.311 },	{ make_pair(KNZ, RNE),	 0.535 },
		{ make_pair(KNZ, SOG),	-0.201 },	{ make_pair(KNZ, HNE),	 0.199 },	{ make_pair(KNZ, YCZ),	-0.348 },
		{ make_pair(KNZ, FCZ),	 0.175 },	{ make_pair(KNZ, LCD),	 0.680 },	{ make_pair(KNZ, CSG),	 1.063 },

		{ make_pair(KCD,   N),	-0.106 },	{ make_pair(KCD,  CA),	-0.066 },	{ make_pair(KCD,   C),	-0.207 },
		{ make_pair(KCD,   O),	-0.323 },	{ make_pair(KCD, GCA),	-0.111 },	{ make_pair(KCD,  CB),	-0.161 },
		{ make_pair(KCD, KNZ),	-0.003 },	{ make_pair(KCD, KCD),	 0.000 },	{ make_pair(KCD, DOD),	-1.405 },
		{ make_pair(KCD, RNH),	 0.240 },	{ make_pair(KCD, NND),	-0.434 },	{ make_pair(KCD, RNE),	 0.465 },
		{ make_pair(KCD, SOG),	-0.324 },	{ make_pair(KCD, HNE),	 0.022 },	{ make_pair(KCD, YCZ),	-0.563 },
		{ make_pair(KCD, FCZ),	-0.123 },	{ make_pair(KCD, LCD),	 0.345 },	{ make_pair(KCD, CSG),	 0.876 },

		{ make_pair(DOD,   N),	 0.126 },	{ make_pair(DOD,  CA),	 0.284 },	{ make_pair(DOD,   C),	 0.151 },
		{ make_pair(DOD,   O),	 0.171 },	{ make_pair(DOD, GCA),	 0.112 },	{ make_pair(DOD,  CB),	 0.105 },
		{ make_pair(DOD, KNZ),	-1.341 },	{ make_pair(DOD, KCD),	-1.405 },	{ make_pair(DOD, DOD),	 0.000 },
		{ make_pair(DOD, RNH),	-1.088 },	{ make_pair(DOD, NND),	-0.253 },	{ make_pair(DOD, RNE),	-1.002 },
		{ make_pair(DOD, SOG),	-0.266 },	{ make_pair(DOD, HNE),	-0.197 },	{ make_pair(DOD, YCZ),	-0.214 },
		{ make_pair(DOD, FCZ),	 0.249 },	{ make_pair(DOD, LCD),	 0.740 },	{ make_pair(DOD, CSG),	 0.853 },

		{ make_pair(RNH,   N),	 0.069 },	{ make_pair(RNH,  CA),	 0.079 },	{ make_pair(RNH,   C),	 0.047 },
		{ make_pair(RNH,   O),	 0.045 },	{ make_pair(RNH, GCA),	 0.027 },	{ make_pair(RNH,  CB),	 0.111 },
		{ make_pair(RNH, KNZ),	 0.331 },	{ make_pair(RNH, KCD),	 0.240 },	{ make_pair(RNH, DOD),	-1.088 },
		{ make_pair(RNH, RNH),	 0.000 },	{ make_pair(RNH, NND),	-0.050 },	{ make_pair(RNH, RNE),	 0.066 },
		{ make_pair(RNH, SOG),	-0.052 },	{ make_pair(RNH, HNE),	 0.304 },	{ make_pair(RNH, YCZ),	-0.181 },
		{ make_pair(RNH, FCZ),	 0.256 },	{ make_pair(RNH, LCD),	 0.728 },	{ make_pair(RNH, CSG),	 1.097 },

		{ make_pair(NND,   N),	 0.117 },	{ make_pair(NND,  CA),	 0.216 },	{ make_pair(NND,   C),	 0.110 },
		{ make_pair(NND,   O),	 0.115 },	{ make_pair(NND, GCA),	 0.131 },	{ make_pair(NND,  CB),	 0.175 },
		{ make_pair(NND, KNZ),	-0.311 },	{ make_pair(NND, KCD),	-0.434 },	{ make_pair(NND, DOD),	-0.253 },
		{ make_pair(NND, RNH),	-0.050 },	{ make_pair(NND, NND),	 0.000 },	{ make_pair(NND, RNE),	-0.022 },
		{ make_pair(NND, SOG),	 0.002 },	{ make_pair(NND, HNE),	 0.330 },	{ make_pair(NND, YCZ),	-0.008 },
		{ make_pair(NND, FCZ),	 0.326 },	{ make_pair(NND, LCD),	 0.675 },	{ make_pair(NND, CSG),	 1.175 },

		{ make_pair(RNE,   N),	 0.203 },	{ make_pair(RNE,  CA),	 0.203 },	{ make_pair(RNE,   C),	 0.116 },
		{ make_pair(RNE,   O),	 0.002 },	{ make_pair(RNE, GCA),	 0.197 },	{ make_pair(RNE,  CB),	 0.145 },
		{ make_pair(RNE, KNZ),	 0.535 },	{ make_pair(RNE, KCD),	 0.465 },	{ make_pair(RNE, DOD),	-1.002 },
		{ make_pair(RNE, RNH),	 0.066 },	{ make_pair(RNE, NND),	-0.022 },	{ make_pair(RNE, RNE),	 0.000 },
		{ make_pair(RNE, SOG),	-0.013 },	{ make_pair(RNE, HNE),	 0.273 },	{ make_pair(RNE, YCZ),	-0.233 },
		{ make_pair(RNE, FCZ),	 0.173 },	{ make_pair(RNE, LCD),	 0.556 },	{ make_pair(RNE, CSG),	 0.964 },

		{ make_pair(SOG,   N),	-0.033 },	{ make_pair(SOG,  CA),	 0.097 },	{ make_pair(SOG,   C),	 0.028 },
		{ make_pair(SOG,   O),	 0.111 },	{ make_pair(SOG, GCA),	 0.064 },	{ make_pair(SOG,  CB),	 0.077 },
		{ make_pair(SOG, KNZ),	-0.201 },	{ make_pair(SOG, KCD),	-0.324 },	{ make_pair(SOG, DOD),	-0.266 },
		{ make_pair(SOG, RNH),	-0.052 },	{ make_pair(SOG, NND),	 0.002 },	{ make_pair(SOG, RNE),	-0.013 },
		{ make_pair(SOG, SOG),	 0.000 },	{ make_pair(SOG, HNE),	 0.190 },	{ make_pair(SOG, YCZ),	 0.179 },
		{ make_pair(SOG, FCZ),	 0.277 },	{ make_pair(SOG, LCD),	 0.503 },	{ make_pair(SOG, CSG),	 0.642 },

		{ make_pair(HNE,   N),	 0.421 },	{ make_pair(HNE,  CA),	 0.406 },	{ make_pair(HNE,   C),	 0.413 },
		{ make_pair(HNE,   O),	 0.399 },	{ make_pair(HNE, GCA),	 0.451 },	{ make_pair(HNE,  CB),	 0.385 },
		{ make_pair(HNE, KNZ),	 0.199 },	{ make_pair(HNE, KCD),	 0.022 },	{ make_pair(HNE, DOD),	-0.197 },
		{ make_pair(HNE, RNH),	 0.304 },	{ make_pair(HNE, NND),	 0.330 },	{ make_pair(HNE, RNE),	 0.273 },
		{ make_pair(HNE, SOG),	 0.190 },	{ make_pair(HNE, HNE),	 0.000 },	{ make_pair(HNE, YCZ),	 0.149 },
		{ make_pair(HNE, FCZ),	 0.362 },	{ make_pair(HNE, LCD),	 0.749 },	{ make_pair(HNE, CSG),	 0.697 },

		{ make_pair(YCZ,   N),	 0.107 },	{ make_pair(YCZ,  CA),	 0.039 },	{ make_pair(YCZ,   C),	 0.094 },
		{ make_pair(YCZ,   O),	 0.092 },	{ make_pair(YCZ, GCA),	 0.004 },	{ make_pair(YCZ,  CB),	-0.046 },
		{ make_pair(YCZ, KNZ),	-0.348 },	{ make_pair(YCZ, KCD),	-0.563 },	{ make_pair(YCZ, DOD),	-0.214 },
		{ make_pair(YCZ, RNH),	-0.181 },	{ make_pair(YCZ, NND),	-0.008 },	{ make_pair(YCZ, RNE),	-0.233 },
		{ make_pair(YCZ, SOG),	 0.179 },	{ make_pair(YCZ, HNE),	 0.149 },	{ make_pair(YCZ, YCZ),	 0.000 },
		{ make_pair(YCZ, FCZ),	 0.023 },	{ make_pair(YCZ, LCD),	 0.130 },	{ make_pair(YCZ, CSG),	 0.618 },

		{ make_pair(FCZ,   N),	 0.207 },	{ make_pair(FCZ,  CA),	 0.175 },	{ make_pair(FCZ,   C),	 0.218 },
		{ make_pair(FCZ,   O),	 0.185 },	{ make_pair(FCZ, GCA),	 0.232 },	{ make_pair(FCZ,  CB),	 0.045 },
		{ make_pair(FCZ, KNZ),	 0.175 },	{ make_pair(FCZ, KCD),	-0.123 },	{ make_pair(FCZ, DOD),	 0.249 },
		{ make_pair(FCZ, RNH),	 0.256 },	{ make_pair(FCZ, NND),	 0.326 },	{ make_pair(FCZ, RNE),	 0.173 },
		{ make_pair(FCZ, SOG),	 0.277 },	{ make_pair(FCZ, HNE),	 0.362 },	{ make_pair(FCZ, YCZ),	 0.023 },
		{ make_pair(FCZ, FCZ),	 0.000 },	{ make_pair(FCZ, LCD),	 0.040 },	{ make_pair(FCZ, CSG),	 0.430 },

		{ make_pair(LCD,   N),	 0.293 },	{ make_pair(LCD,  CA),	 0.235 },	{ make_pair(LCD,   C),	 0.325 },
		{ make_pair(LCD,   O),	 0.294 },	{ make_pair(LCD, GCA),	 0.504 },	{ make_pair(LCD,  CB),	 0.137 },
		{ make_pair(LCD, KNZ),	 0.680 },	{ make_pair(LCD, KCD),	 0.345 },	{ make_pair(LCD, DOD),	 0.740 },
		{ make_pair(LCD, RNH),	 0.728 },	{ make_pair(LCD, NND),	 0.675 },	{ make_pair(LCD, RNE),	 0.556 },
		{ make_pair(LCD, SOG),	 0.503 },	{ make_pair(LCD, HNE),	 0.749 },	{ make_pair(LCD, YCZ),	 0.130 },
		{ make_pair(LCD, FCZ),	 0.040 },	{ make_pair(LCD, LCD),	 0.000 },	{ make_pair(LCD, CSG),	 0.405 },

		{ make_pair(CSG,   N),	 0.173 },	{ make_pair(CSG,  CA),	 0.264 },	{ make_pair(CSG,   C),	 0.190 },
		{ make_pair(CSG,   O),	 0.228 },	{ make_pair(CSG, GCA),	 0.569 },	{ make_pair(CSG,  CB),	 0.306 },
		{ make_pair(CSG, KNZ),	 1.063 },	{ make_pair(CSG, KCD),	 0.876 },	{ make_pair(CSG, DOD),	 0.853 },
		{ make_pair(CSG, RNH),	 1.097 },	{ make_pair(CSG, NND),	 1.175 },	{ make_pair(CSG, RNE),	 0.964 },
		{ make_pair(CSG, SOG),	 0.642 },	{ make_pair(CSG, HNE),	 0.697 },	{ make_pair(CSG, YCZ),	 0.618 },
		{ make_pair(CSG, FCZ),	 0.430 },	{ make_pair(CSG, LCD),	 0.405 },	{ make_pair(CSG, CSG),	 0.000 }
};

/**
 * Given an atom, this method returns the corresponding atomType_t.
 * @param atm	the input atom
 * @return		one of the 18 atom types or OTHER_ATM if the atom type is unknown
 */
atomType_t AtomicContactEnergies::resolveAtomType(atom const & atm) {
	string residue = atm.residue_name;
	string atomName = trim(atm.atom_name);
	switch (str2int(atomName.c_str())) {
	case str2int("C"):
		return C;
		break;
	case str2int("CA"):
		if (residue != "GLY")
			return CA;
		else
			return GCA;
		break;
	case str2int("O"):
		return O;
		break;
	case str2int("N"):
		return N;
		break;
	default:
		switch (str2int(residue.c_str())) {
		case str2int("ALA"):
			if (atomName == "CB")
				return CB;
			else
				return OTHER_ATM;
			break;
		case str2int("ARG"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CG")
				return FCZ;
			else if (atomName == "CD" || atomName == "NE")
				return RNE;
			else if (atomName == "CZ" || atomName == "NH1" || atomName == "NH2")
				return RNH;
			else
				return OTHER_ATM;
			break;
		case str2int("ASN"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CG" || atomName == "ND2" || atomName == "OD1")
				return NND;
			else
				return OTHER_ATM;
			break;
		case str2int("ASP"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CG" || atomName == "OD1" || atomName == "OD2")
				return DOD;
			else
				return OTHER_ATM;
			break;
		case str2int("CYS"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "SG")
				return CSG;
			else
				return OTHER_ATM;
			break;
		case str2int("GLN"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CG")
				return FCZ;
			else if (atomName == "CD" || atomName == "NE2" || atomName == "OE1")
				return NND;
			else
				return OTHER_ATM;
			break;
		case str2int("GLU"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CD" || atomName == "OE1" || atomName == "OE2")
				return DOD;
			else if (atomName == "CG")
				return FCZ;
			else
				return OTHER_ATM;
			break;
		case str2int("GLY"):
			if (atomName == "CA")
				return GCA;
			else
				return OTHER_ATM;
			break;
		case str2int("HIS"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CD2" || atomName == "CE1" || atomName == "CG"
					|| atomName == "ND1" || atomName == "NE2")
				return HNE;
			else
				return OTHER_ATM;
			break;
		case str2int("ILE"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CG1")
				return FCZ;
			else if (atomName == "CD" || atomName == "CG2")
				return LCD;
			else
				return OTHER_ATM;
			break;
		case str2int("LEU"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CG")
				return FCZ;
			else if (atomName == "CD1" || atomName == "CD2")
				return LCD;
			else
				return OTHER_ATM;
			break;
		case str2int("LYS"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CG")
				return FCZ;
			else if (atomName == "CD")
				return KCD;
			else if (atomName == "CE" || atomName == "NZ")
				return KNZ;
			else
				return OTHER_ATM;
			break;
		case str2int("MET"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CG" || atomName == "SD")
				return FCZ;
			else if (atomName == "CE")
				return LCD;
			else
				return OTHER_ATM;
			break;
		case str2int("PHE"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CD1" || atomName == "CD2" || atomName == "CE1"
					|| atomName == "CE2" || atomName == "CG" || atomName == "CZ")
				return FCZ;
			else
				return OTHER_ATM;
			break;
		case str2int("PRO"):
			if (atomName == "CB" || atomName == "CD" || atomName == "CG")
				return CB;
			else
				return OTHER_ATM;
			break;
		case str2int("SER"):
			if (atomName == "CB" || atomName == "OG")
				return SOG;
			else
				return OTHER_ATM;
			break;
		case str2int("THR"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CG2")
				return FCZ;
			else if (atomName == "OG1")
				return SOG;
			else
				return OTHER_ATM;
			break;
		case str2int("TRP"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CD1" || atomName == "CD2"
					|| atomName == "CE2" || atomName == "CE3"
					|| atomName == "CG"	 || atomName == "CH2"
					|| atomName == "CZ2" || atomName == "CZ3")
				return FCZ;
			else if (atomName == "NE1")
				return HNE;
			else
				return OTHER_ATM;
			break;
		case str2int("TYR"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CD1" || atomName == "CD2" || atomName == "CG")
				return FCZ;
			else if (atomName == "OH")
				return SOG;
			else if (atomName == "CE1" || atomName == "CE2" || atomName == "CZ")
				return YCZ;
			else
				return OTHER_ATM;
			break;
		case str2int("VAL"):
			if (atomName == "CB")
				return CB;
			else if (atomName == "CG1" || atomName == "CG2")
				return LCD;
			else
				return OTHER_ATM;
			break;
		default:
			return OTHER_ATM;
		}
	}
};
