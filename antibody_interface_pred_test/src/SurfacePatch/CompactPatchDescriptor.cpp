/*
 * CompactPatchDescriptor.cpp
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */
/**
 * Implementation of the CompatchPatchDescriptor object.
 */
#include "CompactPatchDescriptor.h"
CompactPatchDescriptor::CompactPatchDescriptor() :
	ID(0), /* curvature(-1), */ isInterface(false) { };
CompactPatchDescriptor::CompactPatchDescriptor(CompactPatchDescriptor const & cpd) :
		isInterface(cpd.isInterface), ID(cpd.ID), center(cpd.center),

//		patchNormal(cpd.patchNormal),
//		curvature(cpd.curvature),
//		surfaceI(cpd.surfaceI),
//		potentials_posI(cpd.potentials_posI), potentials_negI(cpd.potentials_negI),
//		hydrophobicity_posI(cpd.hydrophobicity_posI),
//		hydrophobicity_negI(cpd.hydrophobicity_negI),

		BLAM930101_posI(cpd.BLAM930101_posI),
		BLAM930101_negI(cpd.BLAM930101_negI),
		BIOV880101_posI(cpd.BIOV880101_posI),
		BIOV880101_negI(cpd.BIOV880101_negI),
		MAXF760101_I(cpd.MAXF760101_I),
		TSAJ990101_I(cpd.TSAJ990101_I),
		NAKH920108_I(cpd.NAKH920108_I),
		CEDJ970104_I(cpd.CEDJ970104_I),
		LIFS790101_I(cpd.LIFS790101_I),
		MIYS990104_posI(cpd.MIYS990104_posI),
		MIYS990104_negI(cpd.MIYS990104_negI) { };
/**
 * Copy assignment operator
 */
CompactPatchDescriptor & CompactPatchDescriptor::operator=(CompactPatchDescriptor const & cpd) {
	if (this != &cpd) { // protect against invalid self-assignment
		this->isInterface = cpd.isInterface;
		this->ID = cpd.ID;
		this->center = cpd.center;
//		this->surfaceI = cpd.surfaceI;
//		this->potentials_posI = cpd.potentials_posI;
//		this->potentials_negI = cpd.potentials_negI;
//		this->curvature = cpd.curvature;
//		this->patchNormal = cpd.patchNormal;
//		this->hydrophobicity_posI = cpd.hydrophobicity_posI;
//		this->hydrophobicity_negI = cpd.hydrophobicity_negI;

		this->BLAM930101_posI = cpd.BLAM930101_posI;
		this->BLAM930101_negI = cpd.BLAM930101_negI;
		this->BIOV880101_posI = cpd.BIOV880101_posI;
		this->BIOV880101_negI = cpd.BIOV880101_negI;
		this->MAXF760101_I = cpd.MAXF760101_I;
		this->TSAJ990101_I = cpd.TSAJ990101_I;
		this->NAKH920108_I = cpd.NAKH920108_I;
		this->CEDJ970104_I = cpd.CEDJ970104_I;
		this->LIFS790101_I = cpd.LIFS790101_I;
		this->MIYS990104_posI = cpd.MIYS990104_posI;
		this->MIYS990104_negI = cpd.MIYS990104_negI;
	}
	return *this;
}

CompactPatchDescriptor::CompactPatchDescriptor(SurfacePatch const & sp) {
	this->isInterface = sp.interface;
	this->ID = sp.ID;
	this->center = (point3D(sp.patchCenter.ix, sp.patchCenter.iy, sp.patchCenter.iz) / sp.resolution) - sp.translation;
//	this->surfaceI = sp.surfaceD->invariants;
//	this->potentials_posI = sp.potentials_posD->invariants;
//	this->potentials_negI = sp.potentials_negD->invariants;
//	this->curvature = sp.curvature;
//	this->patchNormal = sp.patchNormal;
//	this->hydrophobicity_posI = sp.hydrophobicity_posD->invariants;
//	this->hydrophobicity_negI = sp.hydrophobicity_negD->invariants;

	this->BLAM930101_posI = sp.BLAM930101_posD->invariants;
	this->BLAM930101_negI = sp.BLAM930101_negD->invariants;
	this->BIOV880101_posI = sp.BIOV880101_posD->invariants;
	this->BIOV880101_negI = sp.BIOV880101_negD->invariants;
	this->MAXF760101_I = sp.MAXF760101_D->invariants;
	this->TSAJ990101_I = sp.TSAJ990101_D->invariants;
	this->NAKH920108_I = sp.NAKH920108_D->invariants;
	this->CEDJ970104_I = sp.CEDJ970104_D->invariants;
	this->LIFS790101_I = sp.LIFS790101_D->invariants;
	this->MIYS990104_posI = sp.MIYS990104_posD->invariants;
	this->MIYS990104_negI = sp.MIYS990104_negD->invariants;
}
