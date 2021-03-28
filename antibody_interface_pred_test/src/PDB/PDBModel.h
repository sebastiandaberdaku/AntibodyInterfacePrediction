/**
 * In this header we define the PDBModel class and relative methods.
 */

#ifndef PDBMODEL_H_
#define PDBMODEL_H_

#include "../Geometry/point3D.h"
#include "../Atom/atom.h"
#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <vector>
#include "../utils/hash.h"

typedef enum lineType {
	ATOM, HETATM, OTHER_LN, TER, REMARK, END
} lineType;

class PDBModel {
public:
	size_t line_num;
	std::vector<atom> atomsInModel;
	std::string filename;
	std::string header;
	bool no_hydrogen;
	bool no_hetatm;

	static std::unordered_map<pair<string,string>, float> atom_radii;

	PDBModel(string const & filename, string const & inname_radii, bool no_hydrogen, bool no_hetatm);
	PDBModel(PDBModel const & model);
	/**
	 * Copy assignment operator
	 */
	PDBModel & operator=(PDBModel const & model);

private:
	lineType getLineType(string const & line);
	void loadFile(string const & filename);
	void loadAtomRadii(string const & inname_radii);
	void parseLine(string const & line);
};
#endif /* PDBMODEL_H_ */
