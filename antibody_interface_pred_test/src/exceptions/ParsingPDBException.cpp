/*
 * ParsingPDBException.cpp
 *
 *  Created on: Jan 6, 2014
 *      Author: daberdaku
 */
#include "ParsingPDBException.h"

/**
 * Constructor of the exception.
 * \param txt	Text of the exception.
 * \param mthd	Text identifying the exception's source,
 * 				i.e. the method that threw the exception.
 * \param cs	Text describing the exception's cause.
 *
 */
ParsingPDBException::ParsingPDBException(const std::string& txt =
		"undefined exception", const std::string& mthd = "undefined",
		const std::string& cs = "Generic ParsingPDBException") :
		std::exception(), errorText(txt), methodText(mthd), causeText(cs) {
}

/**
 * Destructor.
 */
ParsingPDBException::~ParsingPDBException() throw () {
}

/**
 * Returns the error text.
 */
std::string ParsingPDBException::error() const {
	return ("Error: " + errorText);
}

/**
 * Returns the method throwing the exception.
 */
std::string ParsingPDBException::method() const {
	if (methodText == "undefined")
		return " ";
	else
		return ("Method: " + methodText);
}

/**
 * Returns the possible cause of the exception.
 */
std::string ParsingPDBException::cause() const {
	return ("Possible cause: " + causeText);
}

/**
 * This method overrides the one implemented in the exception class.
 */
const char* ParsingPDBException::what() const throw () {
	static std::string ex;
	ex = errorText + " -- " + methodText + " -- " + causeText;
	return ex.c_str();
}
