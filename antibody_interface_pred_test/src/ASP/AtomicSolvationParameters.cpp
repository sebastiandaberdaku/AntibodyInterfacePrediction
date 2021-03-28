/*
 * AtomicSolvationParameters.cpp
 *
 *  Created on: 01/mar/2015
 *      Author: sebastian
 */

#include "AtomicSolvationParameters.h"


/**
 *  Zhou H, Zhou Y: Stability scale and atomic solvation parameters extracted from 1023
 *  mutation experiments. Proteins 2002, 49(4):483-492.
 *  Values from the first column of table IV. (unit cal / (mol * Å^2))
 */
const std::map<std::string, float> AtomicSolvationParameters::asp_zhou_1 = {
	{"C",	21}, {"S",	14}, {"O/N",	4}, {"O-",	3}, {"N+",	2}
};
/**
 *  Zhou H, Zhou Y: Stability scale and atomic solvation parameters extracted from 1023
 *  mutation experiments. Proteins 2002, 49(4):483-492.
 *  Values from the second column of table IV. (unit cal / (mol * Å^2))
 */
const std::map<std::string, float> AtomicSolvationParameters::asp_zhou_2 = {
	{"C",	25}, {"S",	18}, {"O/N",	12}, {"O-",	14}, {"N+",	28}
};

/**
 *	Wesson L, Eisenberg D: Atomic solvation parameters applied to molecular
 *	dynamics of proteins in solution. Protein Sci. Feb 1992; 1(2): 227–235.
 *	Values from the first column of table 3. (unit cal / (mol * Å^2))
 */
const std::map<std::string, float> AtomicSolvationParameters::asp_wesson_1 = {
	{ "C", 4 }, { "S", -17 }, { "O/N", -113 }, { "O-", -166 }, { "N+", -169 }
};
/**
 *	Wesson L, Eisenberg D: Atomic solvation parameters applied to molecular
 *	dynamics of proteins in solution. Protein Sci. Feb 1992; 1(2): 227–235.
 *	Values from the second column of table 3. (unit cal / (mol * Å^2))
 */
const std::map<std::string, float> AtomicSolvationParameters::asp_wesson_2 = {
	{ "C", 12 }, { "S", -18 }, { "O/N", -116 }, { "O-", -175 }, { "N+", -186 }
};
/**
 *	Eisenberg D, McLachlan A: Solvation energy in protein folding and binding.
 *	Nature 319, 199 - 203 (16 January 1986)
 *	(unit cal / (mol * Å^2))
 */
const std::map<std::string, float> AtomicSolvationParameters::asp_eisenberg = {
	{ "C", 16 }, { "S", 21 }, { "O/N", -6 }, { "O-", -24 }, { "N+", -50 }
};



