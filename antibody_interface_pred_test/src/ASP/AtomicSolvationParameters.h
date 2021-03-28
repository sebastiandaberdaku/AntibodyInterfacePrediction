/*
 * AtomicSolvationParameters.h
 *
 *  Created on: 21/ott/2014
 *      Author: sebastian
 */
/**
 * In this header we define different Atomic Solvation Parameters scales.
 */
#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_ASP_ATOMICSOLVATIONPARAMETERS_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_ASP_ATOMICSOLVATIONPARAMETERS_H_

#include <map>

typedef struct AtomicSolvationParameters {
	static const std::map<std::string, float> asp_zhou_1;
	static const std::map<std::string, float> asp_zhou_2;
	static const std::map<std::string, float> asp_wesson_1;
	static const std::map<std::string, float> asp_wesson_2;
	static const std::map<std::string, float> asp_eisenberg;
} AtomicSolvationParameters;

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_ASP_ATOMICSOLVATIONPARAMETERS_H_ */
