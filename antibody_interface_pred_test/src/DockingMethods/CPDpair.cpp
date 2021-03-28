/*
 * CPDpair.cpp
 *
 *  Created on: 25/feb/2015
 *      Author: sebastian
 */

#include "CPDpair.h"
#include <fstream>
/**
 * Constructors
 */
CPDpair::CPDpair() :
		cpd1(NULL), cpd2(NULL), score_surface(0), score_electrostatics(0), score_combined(0) { }

CPDpair::CPDpair(CPDpair const & p) :
		cpd1(p.cpd1), cpd2(p.cpd2), score_surface(p.score_surface),
		score_electrostatics(p.score_electrostatics), score_combined(p.score_combined) { }

CPDpair::CPDpair(CompactPatchDescriptor const & cpd1, CompactPatchDescriptor const & cpd2):
		cpd1(&cpd1), cpd2(&cpd2) {

//	double ds = distanceEuclidean(cpd1.surfaceI, cpd2.surfaceI);
//	double de = sqrt(s_distanceEuclidean(cpd1.potentials_posI, cpd2.potentials_negI)
//			+ s_distanceEuclidean(cpd1.potentials_negI, cpd2.potentials_posI));
//
//	score_surface = 1.0 / (1.0 + ds);
//	score_electrostatics = 1.0 / (1.0 + de);
//	score_combined = 0.0;
}

/**
 * Copy assignment operator
 */
CPDpair & CPDpair::operator=(CPDpair const & p) {
	if (this != &p) {
		this->cpd1 = p.cpd1;
		this->cpd2 = p.cpd2;
		this->score_surface = p.score_surface;
		this->score_electrostatics = p.score_electrostatics;
		this->score_combined = p.score_combined;
	}
	return *this;
}
ostream & operator<<(ostream &os, CPDpair const & p) {
//	os << "index 1: " << p.cpd1->ID << " index 2: " << p.cpd2->ID
//			<< " score_surface: " << p.score_surface
//			<< " score_electrostatics: " << p.score_electrostatics << "\n";
//	return os;


	os.precision(std::numeric_limits<double>::digits10 + 2);
	size_t c = 0;
//	for (auto const & d : p.cpd1->surfaceI)
//		os << "\t" << ++c << ":" << d;
//	for (auto const & d : p.cpd1->potentials_posI)
//		os << "\t" << ++c << ":" << d;
//	for (auto const & d : p.cpd1->potentials_negI)
//		os << "\t" << ++c << ":" << d;
//	for (auto const & d : p.cpd1->hydrophobicity_posI)
//		os << "\t" << ++c << ":" << d;
//	for (auto const & d : p.cpd1->hydrophobicity_negI)
//		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->BLAM930101_posI)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->BLAM930101_negI)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->BIOV880101_posI)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->BIOV880101_negI)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->MAXF760101_I)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->TSAJ990101_I)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->NAKH920108_I)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->CEDJ970104_I)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->LIFS790101_I)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->MIYS990104_posI)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd1->MIYS990104_negI)
		os << "\t" << ++c << ":" << d;

//	for (auto const & d : p.cpd2->surfaceI)
//		os << "\t" << ++c << ":" << d;
//	for (auto const & d : p.cpd2->potentials_posI)
//		os << "\t" << ++c << ":" << d;
//	for (auto const & d : p.cpd2->potentials_negI)
//		os << "\t" << ++c << ":" << d;
//	for (auto const & d : p.cpd2->hydrophobicity_posI)
//		os << "\t" << ++c << ":" << d;
//	for (auto const & d : p.cpd2->hydrophobicity_negI)
//		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->BLAM930101_posI)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->BLAM930101_negI)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->BIOV880101_posI)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->BIOV880101_negI)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->MAXF760101_I)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->TSAJ990101_I)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->NAKH920108_I)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->CEDJ970104_I)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->LIFS790101_I)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->MIYS990104_posI)
		os << "\t" << ++c << ":" << d;
	for (auto const & d : p.cpd2->MIYS990104_negI)
		os << "\t" << ++c << ":" << d;

	return os;
}

