/**
 * Implementation of the PDBModel class.
 */

#include "PDBModel.h"

#include "../utils/makeDirectory.h"
#include <parallel/algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "../exceptions/ParsingPDBException.h"

using namespace std;

std::unordered_map<pair<string, string>, float> PDBModel::atom_radii;

/**
 * Constructor
 */
PDBModel::PDBModel(string const & filename, string const & inname_radii, bool no_hydrogen, bool no_hetatm) :
		filename(filename), line_num(0), no_hydrogen(no_hydrogen), no_hetatm(no_hetatm) {
#pragma omp critical
	{
	if (atom_radii.empty())
		loadAtomRadii(inname_radii);
	}
	loadFile(filename);
}
/**
 * Copy constructor
 */
PDBModel::PDBModel(PDBModel const & model) : atomsInModel(model.atomsInModel),
		line_num(model.line_num), filename(model.filename), header(model.header),
		no_hydrogen(model.no_hydrogen), no_hetatm(model.no_hetatm) { }

/**
 * Copy assignment operator
 */
PDBModel & PDBModel::operator=(PDBModel const & model) {
	if (this != &model) {
		this->atomsInModel = model.atomsInModel; //std::vector::operator= (c++11)
		this->line_num = model.line_num;
		this->filename = model.filename;
		this->header = model.header;
	}
	return *this;
}

/**
 * Parses one line from a PDB file. This method gets a string containing the line
 * to be parsed as an input. Each ATOM line record is added to the atomsInModel list.
 *
 * \param line	The input line.
 */
void PDBModel::parseLine(string const & line) {

	lineType lt = getLineType(line);
	if (lt == ATOM || lt == HETATM) {
		/** From http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
		 *
		ATOM

		The ATOM records present the atomic coordinates for standard amino acids and
		nucleotides. They also present the occupancy and temperature factor for each atom.
		Non-polymer chemical coordinates use the HETATM record type. The element symbol is
		always present on each ATOM record; charge is optional.

		Record Format

		COLUMNS        DATA TYPE       CONTENTS
		--------------------------------------------------------------------------------
		 1 -  6        Record name     "ATOM  "
		 7 - 11        Integer         Atom serial number.
		12			   Empty		   Unused
		13 - 16        Atom            Atom name.
		17             Character       Alternate location indicator.
		18 - 20        Residue name    Residue name.
		21			   Empty		   Unused
		22             Character       Chain identifier.
		23 - 26        Integer         Residue sequence number.
		27             AChar           Code for insertion of residues.
		28 - 30		   Empty		   Unused
		31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
		39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
		47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
		55 - 60        Real(6.2)       Occupancy.
		61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
		67 - 72		   Empty		   Unused
		73 - 76        LString(4)      Segment identifier, left-justified.
		77 - 78        LString(2)      Element symbol, right-justified.
		79 - 80        LString(2)      Charge on the atom.

		Example:
		         1         2         3         4         5         6         7         8
		12345678901234567890123456789012345678901234567890123456789012345678901234567890
		ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
		ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C
		ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34      A1   C
		ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65      A1   O
		ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88      A1   C
		ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41      A1   C
		ATOM    151  CG1AVAL A  25      28.870  17.401  57.336  0.28 12.64      A1   C
		ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11      A1   C
		ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58      A1   C
		ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25      A1   C
		ATOM   1660  CE2 TYR B 205      43.549  -4.115  45.779  1.00 19.60           C
		ATOM    547 HD23 LEU A  34      52.930  48.946   9.332  1.00  5.92           H
		 */

		/*
		 * Tokenize record line in the composing entries
		 */
		istringstream atom_record(line);

		char field_name[7] = {'\0'};
		atom_record.get(field_name, 7); 	// 1 -  6        Record name	"ATOM  "
		char atom_number[6] = {'\0'};
		atom_record.get(atom_number, 6);	// 7 - 11        Integer		Atom serial number.
		char empty;
		atom_record.get(empty);				// 12			 Empty			Unused
		char atom_name[5] = {'\0'};
		atom_record.get(atom_name, 5);		// 13 - 16       Atom			Atom name.
		char alt;
		atom_record.get(alt);				// 17            Character		Alternate location indicator.
		char residue_name[4] = {'\0'};
		atom_record.get(residue_name, 4);	// 18 - 20       Residue name	Residue name.
		atom_record.get(empty);				// 21			 Empty			Unused
		char chain_ID;
		atom_record.get(chain_ID);			// 22            Character		Chain identifier.
		char residue_number[5] = {'\0'};
		atom_record.get(residue_number, 5); // 23 - 26       Integer		Residue sequence number.
		char AChar;
		atom_record.get(AChar);				// 27            AChar			Code for insertion of residues.
		char emptyL[4] = {'\0'};
		atom_record.get(emptyL, 4);			// 28 - 30		 Empty			Unused
		char x[9] = {'\0'};
		atom_record.get(x, 9);				// 31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
		char y[9] = {'\0'};
		atom_record.get(y, 9);				// 39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
		char z[9] = {'\0'};
		atom_record.get(z, 9);				// 47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
		char occupancy[7] = {'\0'};
		atom_record.get(occupancy, 7);		// 55 - 60        Real(6.2)       Occupancy.
		char tempFactor[7] = {'\0'};
		atom_record.get(tempFactor, 7);		// 61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
		char emptyLL[7] = {'\0'};
		atom_record.get(emptyLL, 7); 		// 67 - 72		   Empty		   Unused
		char segment_ID[5] = {'\0'};
		atom_record.get(segment_ID, 5);		// 73 - 76        LString(4)      Segment identifier, left-justified.
		char element_symbol[3] = {'\0'};
		atom_record.get(element_symbol, 3);	// 77 - 78        LString(2)      Element symbol, right-justified.
		char charge[3] = {'\0'};
		atom_record.get(charge, 3); 		// 79 - 80        LString(2)      Charge on the atom.

//		if (atom_record.fail()) {
//			ostringstream error;
//			error << "Failed parsing \"" << field_name << "\" line, number: " << line_num
//					<< " in PDB file: " << filename << ".\n";
//			throw ParsingPDBException(error.str(), "PDBModel::parseLine", "Incorrect PDB file format. ");
//		}
		atom atm(field_name, atoi(atom_number), string(atom_name), string(residue_name), chain_ID,
				atoi(residue_number), atof(x), atof(y), atof(z), atof(occupancy), atof(tempFactor));

		if (no_hydrogen && atm.is_Hydrogen())
			return;
		if (no_hetatm && atm.is_HETATM())
			return;

		auto radius = atom_radii.find(make_pair(trim(residue_name), trim(atom_name)));
		if (radius == atom_radii.end() || floatCompare(radius->second, 0)) {
			cout << "Could not assign radius to: " << field_name << " " << atom_number << " "<< atom_name << " " << residue_name<< endl;
			atm.radius = 0;
//			return;
		} else {
			atm.radius = radius->second;
		}

		atomsInModel.push_back(atm);
	} else {

		header += line + "\n";

		if (lt == REMARK) {
			return;

//			istringstream remark_record(line);
//
//			string field_name;
//			string type;
//			string text;
//
//			remark_record >> field_name >> type >> ws;
//			getline(remark_record, text);
//
//		if (type != "1" && type != "6")
//			return;
//		if (text.empty())
//			return;
//		cout << text << endl;
		}

	}
} /* parseLine() */

/**
 * Loads the PDB file. This method parses the input PDB file one line at a time.
 * Each line is parsed by the parseLine() method.
 *
 * \param filename	Input PDB filename.
 */
void PDBModel::loadFile(string const & filename) {
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
		throw ParsingPDBException("Failed opening input PDB file: " + filename, "PDBModel::loadFile", error_message);
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

/** Decides the record type contained in the current line.
 * \param line	The line to be analyzed.
 * \return	The type of line.
 */
lineType PDBModel::getLineType(string const & line) {
	istringstream sline(line);
	string word;
	sline >> word;
//	if (sline.fail()) {
//		stringstream cause;
//		cause << "Unknown line type at line: " << line << "in PDB file: "
//				<< filename << ".";
//		throw ParsingPDBException("\"" + word + "\" is not a known line type.",
//				"PDBModel::getLineType",
//				"Incorrect PDB file format. " + cause.str());
//	}
	if (word.find("ATOM") != string::npos)
		return ATOM;
	else if (word.find("HETATM") != string::npos)
		return HETATM;
	else if (word.find("END") != string::npos)
		return END;
	else if (word.find("REMARK") != string::npos)
		return REMARK;
	else if (word.find("TER") != string::npos)
		return TER;
	else
		return OTHER_LN;
} /* getLineType() */


void PDBModel::loadAtomRadii(string const & inname_radii) {
	std::ifstream file_stream;
	file_stream.open(inname_radii);
	if (!file_stream.is_open()) {
		struct stat st = { 0 };
		stat(inname_radii.c_str(), &st);
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
			error_message = "Unknown error occurred while opening file " + inname_radii + ".";
			break;
		}
		throw ParsingPDBException("Failed opening input atom radius file: " + inname_radii, "PDBModel::loadAtomRadii", error_message);
	}
	//read
	string line;
	size_t line_number = 0;
	while (getline(file_stream, line)) {
		++line_number;
		if (line[0] == '#' || line[0] == '\n')
			continue;
		istringstream record(line);
		string res, atm;
		float rad;
		if (!(record >> res >> atm >> rad))
			throw ParsingPDBException("Failed parsing input atom radius file: " + inname_radii, "PDBModel::loadAtomRadii", "Error at line: " + to_string(line_number));
		atom_radii[make_pair(res, atm)] = rad;
	}
	file_stream.close();
} /* loadFile */
