/*
 * atom.cpp
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

/**
 * Implementation of the atom struct.
 */
#include "../Atom/atom.h"

/**
 * Default constructor.
 */
atom::atom() :
		is_hetero(false), atom_number(0), residue_number(0), chain_ID(0),
		x(min_float), y(min_float), z(min_float), charge(min_float), radius(min_float) { }
/**
 * Constructor.
 */
atom::atom(string const & field_name, int atom_number,
		string const & atom_name, string const & residue_name, char chain_ID,
		int residue_number, float x, float y, float z, float charge, float radius) :
		atom_number(atom_number), atom_name(atom_name),
		residue_name(residue_name), chain_ID(chain_ID), residue_number(residue_number),
		x(x), y(y), z(z), charge(charge), radius(radius) {
		this->is_hetero = (field_name == "HETATM");
}
/**
 * Copy constructor.
 */
atom::atom(atom const & atm) :
		is_hetero(atm.is_hetero), atom_number(atm.atom_number), atom_name(atm.atom_name),
		residue_name(atm.residue_name), chain_ID(atm.chain_ID), residue_number(atm.residue_number),
		x(atm.x), y(atm.y), z(atm.z), charge(atm.charge), radius(atm.radius) { }
/**
 * Point to atom constructor
 * @param atm
 * @return
 */
atom::atom(point3D const & pt) :
		is_hetero(true), atom_number(0), residue_number(0), chain_ID(0),
		x(pt.x), y(pt.y), z(pt.z), charge(min_float), radius(min_float) { }

/**
 * Copy assignment operator
 */
atom & atom::operator=(atom const & atm) {
	if (this != &atm) {
		this->is_hetero = atm.is_hetero;
		this->atom_number = atm.atom_number;
		this->atom_name = atm.atom_name;
		this->residue_name = atm.residue_name;
		this->chain_ID == atm.chain_ID;
		this->residue_number = atm.residue_number;
		this->x = atm.x;
		this->y = atm.y;
		this->z = atm.z;
		this->charge = atm.charge;
		this->radius = atm.radius;
	}
	return *this;
}
/**
 * 'Equal to' operator
 */
bool operator==(atom const & lhs, atom const & rhs) {
	return (lhs.is_hetero == rhs.is_hetero
			&& lhs.atom_number == rhs.atom_number
			&& lhs.atom_name == rhs.atom_name
			&& lhs.residue_name == rhs.residue_name
			&& lhs.chain_ID == rhs.chain_ID
			&& lhs.residue_number == rhs.residue_number
			&& floatCompare(lhs.x, rhs.x)
			&& floatCompare(lhs.y, rhs.y)
			&& floatCompare(lhs.z, rhs.z)
			&& floatCompare(lhs.charge, rhs.charge)
			&& floatCompare(lhs.radius, rhs.radius));
}


