/***********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright (c) 2016, Sebastian Daberdaku (sebastian.daberdaku@dei.unipd.it)
 * DIPARTIMENTO DI INGEGNERIA DELL'INFORMAZIONE
 * Università degli Studi di Padova.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/
/*
 * AAindex.h
 *
 * Amino acid feature set
 *  Created on: Jun 1, 2016
 *      Author: Sebastian Daberdaku
 */

#ifndef AAINDEX_H_
#define AAINDEX_H_

#include <map>
#include <unordered_map>

using namespace std;

/**
 * Alpha helix propensity of position 44 in T4 lysozyme (Blaber et al., 1993)
 */
static const unordered_map<string, float> BLAM930101 = {
		{"ALA",  0.96}, {"ARG",  0.77}, {"ASN",  0.39}, {"ASP",  0.42},
		{"CYS",  0.42}, {"GLN",  0.80}, {"GLU",  0.53}, {"GLY",  0.00},
		{"HIS",  0.57}, {"ILE",  0.84}, {"LEU",  0.92}, {"LYS",  0.73},
		{"MET",  0.86}, {"PHE",  0.59}, {"PRO", -2.50}, {"SER",  0.53},
		{"THR",  0.54}, {"TRP",  0.58}, {"TYR",  0.72}, {"VAL",  0.63}
};
/**
 * Information value for accessibility; average fraction 35% (Biou et al., 1988)
 */
static const unordered_map<string, float> BIOV880101 = {
		{"ALA",  16.0}, {"ARG", -70.0}, {"ASN", -74.0}, {"ASP", -78.0},
		{"CYS", 168.0}, {"GLN", -73.0}, {"GLU",-106.0}, {"GLY", -13.0},
		{"HIS",	 50.0}, {"ILE", 151.0}, {"LEU",	145.0}, {"LYS",-141.0},
		{"MET",	124.0}, {"PHE", 189.0}, {"PRO", -20.0}, {"SER",	-70.0},
		{"THR",	-38.0}, {"TRP", 145.0}, {"TYR",	 53.0}, {"VAL",	123.0}
};
/**
 * Normalized frequency of alpha-helix (Maxfield-Scheraga, 1976)
 */
static const unordered_map<string, float> MAXF760101 = {
		{"ALA",  1.43}, {"ARG",  1.18}, {"ASN",  0.64}, {"ASP",  0.92},
		{"CYS",  0.94}, {"GLN",  1.22}, {"GLU",  1.67}, {"GLY",  0.46},
		{"HIS",	 0.98}, {"ILE",  1.04}, {"LEU",	 1.36}, {"LYS",  1.27},
		{"MET",	 1.53}, {"PHE",  1.19}, {"PRO",  0.49}, {"SER",	 0.70},
		{"THR",	 0.78}, {"TRP",  1.01}, {"TYR",	 0.69}, {"VAL",	 0.98}
};
/**
 * Volumes including the crystallographic waters using the ProtOr (Tsai et al., 1999)
 */
static const unordered_map<string, float> TSAJ990101 = {
		{"ALA",  89.3}, {"ARG", 190.3}, {"ASN", 122.4}, {"ASP", 114.4},
		{"CYS", 102.5}, {"GLN", 146.9}, {"GLU", 138.8}, {"GLY",  63.8},
		{"HIS",	157.5}, {"ILE", 163.0}, {"LEU",	163.1}, {"LYS", 165.1},
		{"MET",	165.8}, {"PHE", 190.8}, {"PRO", 121.6}, {"SER",	 94.2},
		{"THR",	119.6}, {"TRP", 226.4}, {"TYR",	194.6}, {"VAL",	138.2}
};
/**
 * AA composition of MEM of multi-spanning proteins (Nakashima-Nishikawa, 1992)
 */
static const unordered_map<string, float> NAKH920108 = {
		{"ALA",  9.36}, {"ARG",  0.27}, {"ASN",  2.31}, {"ASP",  0.94},
		{"CYS",  2.56}, {"GLN",  1.14}, {"GLU",  0.94}, {"GLY",  6.17},
		{"HIS",	 0.47}, {"ILE", 13.73}, {"LEU",	16.64}, {"LYS",  0.58},
		{"MET",	 3.93}, {"PHE", 10.99}, {"PRO",  1.96}, {"SER",	 5.58},
		{"THR",	 4.68}, {"TRP",  2.20}, {"TYR",	 3.13}, {"VAL",	12.43}
};
/**
 * Composition of amino acids in intracellular proteins (percent) (Cedano et al., 1997)
 */
static const unordered_map<string, float> CEDJ970104 = {
		{"ALA", 7.9}, {"ARG", 4.9}, {"ASN", 4.0}, {"ASP", 5.5},
		{"CYS", 1.9}, {"GLN", 4.4}, {"GLU", 7.1}, {"GLY", 7.1},
		{"HIS",	2.1}, {"ILE", 5.2}, {"LEU", 8.6}, {"LYS", 6.7},
		{"MET",	2.4}, {"PHE", 3.9}, {"PRO", 5.3}, {"SER", 6.6},
		{"THR",	5.3}, {"TRP", 1.2}, {"TYR", 3.1}, {"VAL", 6.8}
};
/**
 * Conformational preference for all beta-strands (Lifson-Sander, 1979)
 */
static const unordered_map<string, float> LIFS790101 = {
		{"ALA",  0.92}, {"ARG",  0.93}, {"ASN",  0.60}, {"ASP",  0.48},
		{"CYS",  1.16}, {"GLN",  0.95}, {"GLU",  0.61}, {"GLY",  0.61},
		{"HIS",	 0.93}, {"ILE",  1.81}, {"LEU",	 1.30}, {"LYS",  0.70},
		{"MET",	 1.19}, {"PHE",  1.25}, {"PRO",  0.40}, {"SER",	 0.82},
		{"THR",	 1.12}, {"TRP",  1.54}, {"TYR",	 1.53}, {"VAL",	 1.81}
};
/**
 * Optimized relative partition energies - method C (Miyazawa-Jernigan, 1999)
 */
static const unordered_map<string, float> MIYS990104 = {
		{"ALA", -0.04}, {"ARG",  0.07}, {"ASN",  0.13}, {"ASP",  0.19},
		{"CYS", -0.38}, {"GLN",  0.14}, {"GLU",  0.23}, {"GLY",  0.09},
		{"HIS",	-0.04}, {"ILE", -0.34}, {"LEU",	-0.37}, {"LYS",  0.33},
		{"MET",	-0.30}, {"PHE", -0.38}, {"PRO",  0.19}, {"SER",	 0.12},
		{"THR",	 0.03}, {"TRP", -0.33}, {"TYR",	-0.29}, {"VAL",	-0.29}
};

/**
 * Method that returns the HQI8 values for a given amino acid.
 * @param aa amino acid name
 * @return	the eight high quality indices for the given amino acid
 */
static vector<float> get_HQI8(string const & aa) {
	vector<float> hqi8(8);
	auto hqi80 = BLAM930101.find(aa);
	if (hqi80 != BLAM930101.end())
		hqi8[0] = hqi80->second;
	else
		hqi8[0] = 0;

	auto hqi81 = BIOV880101.find(aa);
	if (hqi81 != BIOV880101.end())
		hqi8[1] = hqi81->second;
	else
		hqi8[1] = 0;

	auto hqi82 = MAXF760101.find(aa);
	if (hqi82 != MAXF760101.end())
		hqi8[2] = hqi82->second;
	else
		hqi8[2] = 0;

	auto hqi83 = TSAJ990101.find(aa);
	if (hqi83 != TSAJ990101.end())
		hqi8[3] = hqi83->second;
	else
		hqi8[3] = 0;

	auto hqi84 = NAKH920108.find(aa);
	if (hqi84 != NAKH920108.end())
		hqi8[4] = hqi84->second;
	else
		hqi8[4] = 0;

	auto hqi85 = CEDJ970104.find(aa);
	if (hqi85 != CEDJ970104.end())
		hqi8[5] = hqi85->second;
	else
		hqi8[5] = 0;

	auto hqi86 = LIFS790101.find(aa);
	if (hqi86 != LIFS790101.end())
		hqi8[6] = hqi86->second;
	else
		hqi8[6] = 0;

	auto hqi87 = MIYS990104.find(aa);
	if (hqi87 != MIYS990104.end())
		hqi8[7] = hqi87->second;
	else
		hqi8[7] = 0;

	return hqi8;
}
#endif /* AAINDEX_H_ */
