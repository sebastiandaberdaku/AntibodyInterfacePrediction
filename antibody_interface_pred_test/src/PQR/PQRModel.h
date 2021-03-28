/**
 * In this header we define the PQRModel class and relative methods.
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_PQR_PQRMODEL_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_PQR_PQRMODEL_H_

#include "../Geometry/point3D.h"
#include "../DockingMethods/DockingPose.h"
#include <vector>
#include "../PDB/PDBModel.h"

class PQRModel {
public:
	size_t line_num;
	std::vector<atom> atomsInModel;
	std::string filename;
	std::string header;
	bool no_hydrogen;
	bool no_hetatm;

	PQRModel(string const & filename, bool no_hydrogen, bool no_hetatm);
	PQRModel(PQRModel const & model);
	/**
	 * Copy assignment operator
	 */
	PQRModel & operator=(PQRModel const & model);
	void merge(PQRModel const & model);
	void transform(DockingPose const & pose);
	void outputPQRFile(string const & filename);
	inline friend std::ostream & operator<<(std::ostream & os, PQRModel const & m) {
		for (auto const & atm : m.atomsInModel)
			os << atm << endl;
		return os;
	}


private:
	lineType getLineType(string const & line);
	void loadFile(string const & filename);
	void parseLine(string const & line);
};
#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_PQR_PQRMODEL_H_ */
