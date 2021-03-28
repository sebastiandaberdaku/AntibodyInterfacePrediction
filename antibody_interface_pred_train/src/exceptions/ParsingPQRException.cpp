/*
 * ParsingPQRException.cpp
 *
 *  Created on: Jan 6, 2014
 *      Author: daberdaku
 */
#include "ParsingPQRException.h"

/**
 * Constructor of the exception.
 * \param txt	Text of the exception.
 * \param mthd	Text identifying the exception's source,
 * 				i.e. the method that threw the exception.
 * \param cs	Text describing the exception's cause.
 *
 */
ParsingPQRException::ParsingPQRException(const std::string& txt =
		"undefined exception", const std::string& mthd = "undefined",
		const std::string& cs = "Generic ParsingPQRException") :
		std::exception(), errorText(txt), methodText(mthd), causeText(cs) {
}

/**
 * Destructor.
 */
ParsingPQRException::~ParsingPQRException() throw () {
}

/**
 * Returns the error text.
 */
std::string ParsingPQRException::error() const {
	return ("Error: " + errorText);
}

/**
 * Returns the method throwing the exception.
 */
std::string ParsingPQRException::method() const {
	if (methodText == "undefined")
		return " ";
	else
		return ("Method: " + methodText);
}

/**
 * Returns the possible cause of the exception.
 */
std::string ParsingPQRException::cause() const {
	return ("Possible cause: " + causeText);
}

/**
 * This method overrides the one implemented in the exception class.
 */
const char* ParsingPQRException::what() const throw () {
	static std::string ex;
	ex = errorText + " -- " + methodText + " -- " + causeText;
	return ex.c_str();
}
