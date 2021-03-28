/*
 * ParsingPDBException.h
 *
 *  Created on: Jan 6, 2014
 *      Author: daberdaku
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_EXCEPTIONS_PARSINGPQREXCEPTION_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_EXCEPTIONS_PARSINGPQREXCEPTION_H_

#include <exception>
#include <iostream>

class ParsingPQRException: public std::exception {
public:
	ParsingPQRException(const std::string& txt, const std::string& mthd, const std::string& cs);
	virtual ~ParsingPQRException() throw();
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
#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_EXCEPTIONS_PARSINGPQREXCEPTION_H_ */
