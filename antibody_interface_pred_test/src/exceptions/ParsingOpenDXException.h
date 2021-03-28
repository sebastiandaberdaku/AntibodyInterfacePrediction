/*
 * ParsingOpenDXException.h
 *
 *  Created on: 11/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_EXCEPTIONS_PARSINGOPENDXEXCEPTION_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_EXCEPTIONS_PARSINGOPENDXEXCEPTION_H_

#include <exception>
#include <iostream>

class ParsingOpenDXException: public std::exception {
public:
	ParsingOpenDXException(const std::string& txt, const std::string& mthd, const std::string& cs);
	virtual ~ParsingOpenDXException() throw();
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


#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_EXCEPTIONS_PARSINGOPENDXEXCEPTION_H_ */
