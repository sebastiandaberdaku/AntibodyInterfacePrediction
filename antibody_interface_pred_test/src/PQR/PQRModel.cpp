/**
 * Implementation of the PQRModel class.
 */

#include "../exceptions/ParsingPQRException.h"
#include "../utils/makeDirectory.h"
#include "../utils/stoc.h"
#include "PQRModel.h"
#include <parallel/algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;
/**
 * Constructor
 */
PQRModel::PQRModel(string const & filename, bool no_hydrogen, bool no_hetatm) :
		filename(filename), line_num(0), no_hydrogen(no_hydrogen), no_hetatm(no_hetatm){
	loadFile(filename);
}
/**
 * Copy constructor
 */
PQRModel::PQRModel(PQRModel const & model) : atomsInModel(model.atomsInModel),
		line_num(model.line_num), filename(model.filename), header(model.header),
		no_hydrogen(model.no_hydrogen), no_hetatm(model.no_hetatm) { }

/**
 * Copy assignment operator
 */
PQRModel & PQRModel::operator=(PQRModel const & model) {
	if (this != &model) {
		this->atomsInModel = model.atomsInModel; //std::vector::operator= (c++11)
		this->line_num = model.line_num;
		this->filename = model.filename;
		this->header = model.header;
	}
	return *this;
}

/**
 * Parses one line from a PQR file. This method gets a string containing the line
 * to be parsed as an input. Each ATOM line record is added to the atomsInModel list.
 *
 * \param line	The input line.
 */
void PQRModel::parseLine(string const & line) {

	lineType lt = getLineType(line);
	if (lt == ATOM || lt == HETATM) {
		/*
		 * Tokenize record line in the composing entries
		 */
		istringstream atom_record(line);

		string field_name, atom_number, residue_name, chain_ID,
				residue_number, x, y, z, charge, radius;
		char atom_name[5] = {'\0'};

		atom_record >> field_name >> atom_number;

		atom_record.get();
		atom_record.get(atom_name, 5);

		atom_record >> residue_name >> ws;

		char ID;
		atom_record.get(ID);
		chain_ID.push_back(ID);

		atom_record >> residue_number >> x >> y >> z >> charge >> radius;

		if (atom_record.fail()) {
			ostringstream error;
			error << "Failed parsing \"" << field_name << "\" line, number: " << line_num
					<< " in PQR file: " << filename << ".\n";
			ostringstream cause;
			cause << "The input PQR file is probably missing the chain ID field, or there are no whitespaces between atom name and residue name, between x and y, or between y and z.\n";
			cause << "Try adding the \"--chain\" option to PDB2PQR to keep the chain id in the output PQR file.\n";
			cause << "Try adding the \"--whitespace\" option to PDB2PQR to insert whitespaces between atom name and residue name, between x and y, and between y and z.\n";
			cause << "If using the PDB2PQR Server, please set the \"Add/keep chain IDs in the PQR file\" and \"Insert whitespaces between atom name and residue name, between x and y, and between y and z\" options in the \"Available options\" section.";
			throw ParsingPQRException(error.str(), "PQRModel::parseLine",
					"Incorrect PQR file format. " + cause.str());
		}
		atom atm(field_name, stoi(atom_number), string(atom_name), residue_name, stoc(chain_ID),
						stoi(residue_number), stof(x), stof(y), stof(z), stof(charge), stof(radius));
		if (no_hydrogen && atm.is_Hydrogen())
			return;
		if (no_hetatm && atm.is_HETATM())
			return;
		atomsInModel.push_back(atm);
	} else {

		header += line + "\n";

		if (lt == REMARK) {

			// REMARK   6 Total charge on this protein: 8.0000 e
			istringstream remark_record(line);

			string field_name;
			string type;
			string text;

			remark_record >> field_name >> type >> ws;
			getline(remark_record, text);

		if (type != "1" && type != "6")
			return;
			if (text.empty())
				return;
			cout << text << endl;
		}

	}
} /* parseLine() */

/**
 * Loads the PQR file. This method parses the input PQR file one line at a time.
 * Each line is parsed by the parseLine() method.
 *
 * \param filename	Input PQR filename.
 */
void PQRModel::loadFile(string const & filename) {
	std::ifstream file_stream;
	file_stream.open(filename);
	if (!file_stream.is_open()) {
		struct stat st = { 0 };
		stat(filename.c_str(), &st);
		string error_message;
		switch (errno) {
		case EACCES:
			error_message = "Search permission is denied for one of the directories in the path prefix of path.";
			break;
		case EBADF:
			error_message = "File descriptor is bad.";
			break;
		case EFAULT:
			error_message = "Pathname points outside your accessible address space.";
			break;
		case ELOOP:
			error_message = "Too many symbolic links were encountered in resolving pathname.";
			break;
		case ENAMETOOLONG:
			error_message = "Path is too long.";
			break;
		case ENOENT:
			error_message = "The file or a directory component in pathname does not exist or is a dangling symbolic link.";
			break;
		case ENOMEM:
			error_message = "Insufficient kernel memory was available.";
			break;
		case ENOTDIR:
			error_message = "A component used as a directory in pathname is not, in fact, a directory.";
			break;
		case EOVERFLOW:
			error_message = "Path refers to a file whose size, inode number, or number of blocks cannot be represented in, respectively, the types off_t, ino_t, or blkcnt_t. This error can occur when, for example, an application compiled on a 32-bit platform without -D_FILE_OFFSET_BITS=64 calls stat() on a file whose size exceeds (1<<31)-1 bytes.";
			break;
		default:
			error_message = "Unknown error occurred while opening file " + filename + ".";
			break;
		}
		throw ParsingPQRException("Failed opening input PQR file: " + filename, "PQRModel::loadFile", error_message);
	}
	//read
	line_num = 0;
	string line;
	while (getline(file_stream, line)) {
		++line_num;
		parseLine(line);
	}
	file_stream.close();
} /* loadFile */

/**
 * Merges the current PQRModel with model.
 */
void PQRModel::merge(PQRModel const & model){
	if(model.atomsInModel.empty())
		return;
	this->atomsInModel.insert(this->atomsInModel.end(), model.atomsInModel.begin(),
			model.atomsInModel.end());
	this->line_num += model.line_num;
}

/**
* Applies a rotation and a translation to the all the atoms in the input vector.
* @param p		the DockingPose to be applied (rotation + translation)
*/
void PQRModel::transform(DockingPose const & p) {
	if (atomsInModel.empty())
		return;
	for (auto & atm : atomsInModel)
		atm.transform(p.rotation, p.translation);
}

void PQRModel::outputPQRFile(string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	string fname = "./output/" + filename + ".pqr";
	file_stream.open(fname);
	if (!file_stream.is_open()) {
		struct stat st = { 0 };
		stat(fname.c_str(), &st);
		string error_message;
		switch (errno) {
		case EACCES:
			error_message = "Search permission is denied for one of the directories in the path prefix of path.";
			break;
		case EBADF:
			error_message = "File descriptor is bad.";
			break;
		case EFAULT:
			error_message = "Pathname points outside your accessible address space.";
			break;
		case ELOOP:
			error_message = "Too many symbolic links were encountered in resolving pathname.";
			break;
		case ENAMETOOLONG:
			error_message = "Path is too long.";
			break;
		case ENOENT:
			error_message = "A directory component in pathname does not exist or is a dangling symbolic link.";
			break;
		case ENOMEM:
			error_message = "Insufficient kernel memory was available.";
			break;
		case ENOTDIR:
			error_message = "A component used as a directory in pathname is not, in fact, a directory.";
			break;
		case EOVERFLOW:
			error_message = "Path refers to a file whose size, inode number, or number of blocks cannot be represented in, respectively, the types off_t, ino_t, or blkcnt_t. This error can occur when, for example, an application compiled on a 32-bit platform without -D_FILE_OFFSET_BITS=64 calls stat() on a file whose size exceeds (1<<31)-1 bytes.";
			break;
		default:
			error_message = "Unknown error occurred while opening file " + fname + ".";
			break;
		}
		throw ParsingPQRException("Failed opening output PQR file: " + fname, "PQRModel::outputPQRFile", error_message);
	}

	file_stream << this->header;
	for (auto const & atm : atomsInModel) {
		file_stream << atm << "\n";
	}

	file_stream.close();
}

/** Decides the record type contained in the current line.
 * \param line	The line to be analyzed.
 * \return	The type of line.
 */
lineType PQRModel::getLineType(string const & line) {
	istringstream sline(line);
	string word;
	sline >> word;
	if (sline.fail()) {
		stringstream cause;
		cause << "Unknown line type at line: " << line << "in PQR file: "
				<< filename << ".";
		throw ParsingPQRException("\"" + word + "\" is not a known line type.",
				"PQRModel::getLineType",
				"Incorrect PQR file format. " + cause.str());
	}

	if (word == "ATOM")
		return ATOM;
	else if (word == "HETATM")
		return HETATM;
	else if (word == "END")
		return END;
	else if (word == "REMARK")
		return REMARK;
	else if (word == "TER")
		return TER;
	else
		return OTHER_LN;
} /* getLineType() */
