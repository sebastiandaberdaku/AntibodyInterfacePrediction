/*
 * ParsingPDBException.h
 *
 *  Created on: Jan 6, 2014
 *      Author: daberdaku
 */

#ifndef PARSING_PDB_EXCEPTION_H_
#define PARSING_PDB_EXCEPTION_H_

#include <exception>
#include <iostream>

class ParsingPDBException: public std::exception {
public:
	ParsingPDBException(const std::string& txt, const std::string& mthd, const std::string& cs);
	virtual ~ParsingPDBException() throw();
	std::string error() const;
	std::string method() const;
	std::string cause() const;
	const char* what() const throw ();

private:
	/**
	 * The text of the exception message.
	 */
	std::string errorText;
	/**
	 * The argument related to this exception.
	 */
	std::string methodText;
	/**
	 * Describes the type of the exception.
	 */
	std::string causeText;
};
#endif /* PARSING_PDB_EXCEPTION_H_ */
