/*
 * ParsingOpenDXException.cpp
 *
 *  Created on: 11/ago/2014
 *      Author: sebastian
 */

#include "ParsingOpenDXException.h"

/**
 * Constructor of the exception.
 * \param txt	Text of the exception.
 * \param mthd	Text identifying the exception's source,
 * 				i.e. the method that threw the exception.
 * \param cs	Text describing the exception's cause.
 *
 */
ParsingOpenDXException::ParsingOpenDXException(const std::string& txt =
		"undefined exception", const std::string& mthd = "undefined",
		const std::string& cs = "Generic ParsingOpenDXException") :
		std::exception(), errorText(txt), methodText(mthd), causeText(cs) {
}

/**
 * Destructor.
 */
ParsingOpenDXException::~ParsingOpenDXException() throw () {
}

/**
 * Returns the error text.
 */
std::string ParsingOpenDXException::error() const {
	return ("Error: " + errorText);
}

/**
 * Returns the method throwing the exception.
 */
std::string ParsingOpenDXException::method() const {
	if (methodText == "undefined")
		return " ";
	else
		return ("Method: " + methodText);
}

/**
 * Returns the possible cause of the exception.
 */
std::string ParsingOpenDXException::cause() const {
	return ("Possible cause: " + causeText);
}

/**
 * This method overrides the one implemented in the exception class.
 */
const char* ParsingOpenDXException::what() const throw () {
	static std::string ex;
	ex = errorText + " -- " + methodText + " -- " + causeText;
	return ex.c_str();
}
