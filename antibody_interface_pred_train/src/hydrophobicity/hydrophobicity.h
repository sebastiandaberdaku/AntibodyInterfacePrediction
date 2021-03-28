/*
 * hydrophobicity.h
 *
 *  Created on: 11/ago/2014
 *      Author: sebastian
 */

/**
 * Hydrophobicity scales are values that define relative hydrophobicity of
 * amino acid residues. The more positive the value, the more hydrophobic
 * are the amino acids located in that region of the protein.
 *
 * A number of different hydrophobicity scales have been developed. There
 * is a clear difference between these scales due to the different methods
 * used to measure hydrophobicity.
 */
#ifndef HYDROPHOBICITY_H_
#define HYDROPHOBICITY_H_

#include <map>
#include <unordered_map>
#include <boost/functional/hash_fwd.hpp>
#include "../utils/hash.h"

using namespace std;
/**
 * Mean fractional area loss (f) [average area buried/standard state area].
 * The hydrophobicity scale by Rose et al. is correlated to the average area
 * of buried amino acids in globular proteins [Rose et al., 1985]. This results
 * in a scale which is not showing the helices of a protein, but rather the
 * surface accessibility.
 * Author(s): Rose G.D., Geselowitz A.R., Lesser G.J., Lee R.H., Zehfus M.H.
 * Reference: Science 229:834-838(1985).
 */
static const std::map<std::string, float> hphob_rose = {
	{"ALA",	0.74}, {"ARG",	0.64}, {"ASN",	0.63}, {"ASP",	0.62},
	{"CYS",	0.91}, {"GLN",	0.62}, {"GLU",	0.62}, {"GLY",	0.72},
	{"HIS",	0.78}, {"ILE",	0.88}, {"LEU",	0.85}, {"LYS",	0.52},
	{"MET",	0.85}, {"PHE",	0.88}, {"PRO",	0.64}, {"SER",	0.66},
	{"THR",	0.70}, {"TRP",	0.85}, {"TYR",	0.76}, {"VAL",	0.86}
};

/**
 * Optimized matching hydrophobicity (OMH).
 * Author(s): Sweet R.M., Eisenberg D.
 * Reference: J. Mol. Biol. 171:479-488(1983).
 */
static const std::map<std::string, float> hphob_sweet = {
	{"ALA",	-0.40}, {"ARG",	-0.59}, {"ASN",	-0.92}, {"ASP",	-1.31},
	{"CYS",	 0.17}, {"GLN",	-0.91}, {"GLU",	-1.22}, {"GLY",	-0.67},
	{"HIS",	-0.64}, {"ILE",	 1.25}, {"LEU",	 1.22}, {"LYS",	-0.67},
	{"MET",	 1.02}, {"PHE",	 1.92}, {"PRO",	-0.49}, {"SER",	-0.55},
	{"THR",	-0.28}, {"TRP",	 0.50}, {"TYR",	 1.67}, {"VAL",	 0.91}
};

/**
 * Hydropathicity.
 * The Kyte-Doolittle scale is widely used for detecting hydrophobic
 * regions in proteins. Regions with a positive value are hydrophobic.
 * This scale can be used for identifying both surface-exposed regions
 * as well as transmembrane regions, depending on the window size used.
 * Short window sizes of 5-7 generally work well for predicting putative
 * surface-exposed regions. Large window sizes of 19-21 are well suited
 * for finding transmembrane domains if the values calculated are above
 * 1.6. These values should be used as a rule of thumb and deviations
 * from the rule may occur.
 * Author(s): Kyte J., Doolittle R.F.
 * Reference: J. Mol. Biol. 157:105-132(1982).
 */
static const std::map<std::string, float> hphob_kyte = {
	{"ALA",	 1.80}, {"ARG",	-4.50}, {"ASN",	-3.50}, {"ASP",	-3.50},
	{"CYS",	 2.50}, {"GLN",	-3.50}, {"GLU",	-3.50}, {"GLY",	-0.40},
	{"HIS",	-3.20}, {"ILE",	 4.50}, {"LEU",	 3.80}, {"LYS",	-3.90},
	{"MET",	 1.90}, {"PHE",	 2.80}, {"PRO",	-1.60}, {"SER",	-0.80},
	{"THR",	-0.70}, {"TRP",	-0.90}, {"TYR",	-1.30}, {"VAL",	 4.20}
};

/**
 * Hydrophobicity (delta G1/2 cal)
 * Author(s): Abraham D.J., Leo A.J.
 * Reference: Proteins: Structure, Function and Genetics 2:130-152(1987).
 */
static const std::map<std::string, float> hphob_abraham = {
	{"ALA",	 0.44}, {"ARG",	-2.42}, {"ASN",	-1.32}, {"ASP",	-0.31},
	{"CYS",	 0.58}, {"GLN",	-0.71}, {"GLU",	-0.34}, {"GLY",	 0.00},
	{"HIS",	-0.01}, {"ILE",	 2.46}, {"LEU",	 2.46}, {"LYS",	-2.45},
	{"MET",	 1.10}, {"PHE",	 2.54}, {"PRO",	 1.29}, {"SER",	-0.84},
	{"THR",	-0.41}, {"TRP",	 2.56}, {"TYR",	 1.63}, {"VAL",	 1.73}
};

/**
 * Hydrophobicity (free energy of transfer to surface in kcal/mole).
 * Author(s): Bull H.B., Breese K.
 * Reference: Arch. Biochem. Biophys. 161:665-670(1974).
 */
static const std::map<std::string, float> hphob_bull = {
	{"ALA",  0.61}, {"ARG",  0.69}, {"ASN",  0.89}, {"ASP",  0.61},
	{"CYS",  0.36}, {"GLN",  0.97}, {"GLU",  0.51}, {"GLY",  0.81},
	{"HIS",  0.69}, {"ILE", -1.45}, {"LEU", -1.65}, {"LYS",  0.46},
	{"MET", -0.66}, {"PHE", -1.52}, {"PRO", -0.17}, {"SER",  0.42},
	{"THR",  0.29}, {"TRP", -1.20}, {"TYR", -1.43}, {"VAL", -0.75}
};

/**
 * Hydrophobicity scale based on free energy of transfer (kcal/mole).
 * Author(s): Guy H.R.
 * Reference: Biophys J. 47:61-70(1985).
 */
static const std::map<std::string, float> hphob_guy = {
	{"ALA",  0.10}, {"ARG",  1.91}, {"ASN",  0.48}, {"ASP",  0.78},
	{"CYS", -1.42},	{"GLN",  0.95},	{"GLU",  0.83},	{"GLY",  0.33},
	{"HIS", -0.50},	{"ILE", -1.13},	{"LEU", -1.18},	{"LYS",  1.40},
	{"MET", -1.59},	{"PHE", -2.12},	{"PRO",  0.73},	{"SER",  0.52},
	{"THR",  0.07}, {"TRP", -0.51}, {"TYR", -0.21}, {"VAL", -1.27}
};

/**
 * Hydrophobicity scale (contact energy derived from 3D data).
 * Author(s): Miyazawa S., Jernigen R.L.
 * Reference: Macromolecules 18:534-552(1985).
*/
static const std::map<std::string, float> hphob_miyazawa = {
	{"ALA",  5.33}, {"ARG",  4.18}, {"ASN",  3.71}, {"ASP",  3.59},
	{"CYS",  7.93}, {"GLN",  3.87}, {"GLU",  3.65}, {"GLY",  4.48},
	{"HIS",  5.10}, {"ILE",  8.83}, {"LEU",  8.47}, {"LYS",  2.95},
	{"MET",  8.95}, {"PHE",  9.03}, {"PRO",  3.87}, {"SER",  4.09},
	{"THR",  4.49}, {"TRP",  7.66}, {"TYR",  5.89}, {"VAL",  7.63}
};

/**
 * Hydrophobicity scale (pi-r).
 * Author(s): Roseman M.A.
 * Reference: J. Mol. Biol. 200:513-522(1988).
 */
static const std::map<std::string, float> hphob_roseman = {
	{"ALA",  0.39}, {"ARG", -3.95}, {"ASN", -1.91}, {"ASP", -3.81},
	{"CYS",  0.25}, {"GLN", -1.30}, {"GLU", -2.91}, {"GLY",  0.00},
	{"HIS", -0.64}, {"ILE",  1.82}, {"LEU",  1.82}, {"LYS", -2.77},
	{"MET",  0.96}, {"PHE",  2.27}, {"PRO",  0.99}, {"SER", -1.24},
	{"THR", -1.00}, {"TRP",  2.13}, {"TYR",  1.47}, {"VAL",  1.30}
};

/**
 * Hydration potential (kcal/mole) at 25øC.
 * Author(s): Wolfenden R.V., Andersson L., Cullis P.M., Southgate C.C.F.
 * Reference: Biochemistry 20:849-855(1981).
 */
static const std::map<std::string, float> hphob_wolfenden = {
	{"ALA",  1.94}, {"ARG", -19.92}, {"ASN", -9.68}, {"ASP", -10.95},
	{"CYS", -1.24}, {"GLN", -9.38}, {"GLU", -10.20}, {"GLY",  2.39},
	{"HIS", -10.27}, {"ILE",  2.15}, {"LEU",  2.28}, {"LYS", -9.52},
	{"MET", -1.48}, {"PHE", -0.76}, {"PRO",  0.00}, {"SER", -5.06},
	{"THR", -4.88}, {"TRP", -5.88}, {"TYR", -6.11}, {"VAL",  1.99}
};

/**
 * Hydrophobic constants derived from HPLC peptide retention times.
 * Author(s): Wilson K.J., Honegger A., Stotzel R.P., Hughes G.J.
 * Reference: Biochem. J. 199:31-41(1981).
 */
static const std::map<std::string, float> hphob_wilson = {
	{"ALA", -0.30}, {"ARG", -1.10}, {"ASN", -0.20}, {"ASP", -1.40},
	{"CYS",  6.30}, {"GLN", -0.20}, {"GLU",  0.00}, {"GLY",  1.20},
	{"HIS", -1.30}, {"ILE",  4.30}, {"LEU",  6.60}, {"LYS", -3.60},
	{"MET",  2.50}, {"PHE",  7.50}, {"PRO",  2.20}, {"SER", -0.60},
	{"THR", -2.20}, {"TRP",  7.90}, {"TYR",  7.10}, {"VAL",  5.90}
};

/**
 * Hydrophobicity indices at ph 3.4 determined by HPLC.
 * Author(s): Cowan R., Whittaker R.G.
 * Reference: Peptide Research 3:75-80(1990).
 */
static const std::map<std::string, float> hphob_cowan = {
	{"ALA",  0.42}, {"ARG", -1.56}, {"ASN", -1.03}, {"ASP", -0.51},
	{"CYS",  0.84}, {"GLN", -0.96}, {"GLU", -0.37}, {"GLY",  0.00},
	{"HIS", -2.28}, {"ILE",  1.81}, {"LEU",  1.80}, {"LYS", -2.03},
	{"MET",  1.18}, {"PHE",  1.74}, {"PRO",  0.86}, {"SER", -0.64},
	{"THR", -0.26}, {"TRP",  1.46}, {"TYR",  0.51}, {"VAL",  1.34}
};

/**
 * Mobilities of amino acids on chromatography paper (RF).
 * Author(s): Aboderin A.A.
 * Reference: Int. J. Biochem. 2:537-544(1971).
 */
static const std::map<std::string, float> hphob_aboderin = {
	{"ALA",  5.10}, {"ARG",  2.00}, {"ASN",  0.60}, {"ASP",  0.70},
	{"CYS",  0.00}, {"GLN",  1.40}, {"GLU",  1.80}, {"GLY",  4.10},
	{"HIS",  1.60}, {"ILE",  9.30}, {"LEU", 10.00}, {"LYS",  1.30},
	{"MET",  8.70}, {"PHE",  9.60}, {"PRO",  4.90}, {"SER",  3.10},
	{"THR",  3.50}, {"TRP",  9.20}, {"TYR",  8.00}, {"VAL",  8.50}
};

/**
 * Proportion of residues 95% buried (in 12 proteins).
 * Author(s): Chothia C.
 * Reference: J. Mol. Biol. 105:1-14(1976).
 */
static const std::map<std::string, float> hphob_chothia = {
	{"ALA",  0.38}, {"ARG",  0.01}, {"ASN",  0.12}, {"ASP",  0.15},
	{"CYS",  0.50}, {"GLN",  0.07}, {"GLU",  0.18}, {"GLY",  0.36},
	{"HIS",  0.17}, {"ILE",  0.60}, {"LEU",  0.45}, {"LYS",  0.03},
	{"MET",  0.40}, {"PHE",  0.50}, {"PRO",  0.18}, {"SER",  0.22},
	{"THR",  0.23}, {"TRP",  0.27}, {"TYR",  0.15}, {"VAL",  0.54}
};

/**
 * The Eisenberg scale is a normalized consensus hydrophobicity scale
 * which shares many features with the other hydrophobicity scales.
 * Author(s): Eisenberg D., Schwarz E., Komarony M., Wall R.
 * Reference: J. Mol. Biol. 179:125-142(1984).
 */
static const std::map<std::string, float> hphob_eisenberg = {
	{"ALA",  0.62}, {"ARG", -2.53}, {"ASN", -0.78}, {"ASP", -0.90},
	{"CYS",  0.29}, {"GLN", -0.85}, {"GLU", -0.74}, {"GLY",  0.48},
	{"HIS", -0.40}, {"ILE",  1.38}, {"LEU",  1.06}, {"LYS", -1.50},
	{"MET",  0.64}, {"PHE",  1.19}, {"PRO",  0.12}, {"SER", -0.18},
	{"THR", -0.05}, {"TRP",  0.81}, {"TYR",  0.26}, {"VAL",  1.08}
};

/**
 * Hopp and Woods developed their hydrophobicity scale for identification
 * of potentially antigenic sites in proteins. This scale is basically a
 * hydrophilic index where apolar residues have been assigned negative
 * values. Antigenic sites are likely to be predicted when using a window
 * size of 7
 * Author(s): Hopp T.P., Woods K.R.
 * Reference: Proc. Natl. Acad. Sci. U.S.A. 78:3824-3828(1981).
 */
static const std::map<std::string, float> hphob_hopp = {
	{"ALA", -0.50}, {"ARG",  3.00}, {"ASN",  0.20}, {"ASP",  3.00},
	{"CYS", -1.00}, {"GLN",  0.20}, {"GLU",  3.00}, {"GLY",  0.00},
	{"HIS", -0.50}, {"ILE", -1.80}, {"LEU", -1.80}, {"LYS",  3.00},
	{"MET", -1.30}, {"PHE", -2.50}, {"PRO",  0.00}, {"SER",  0.30},
	{"THR", -0.40}, {"TRP", -3.40}, {"TYR", -2.30}, {"VAL", -1.50}
};

/**
 * Average surrounding hydrophobicity.
 * Author(s): Manavalan P., Ponnuswamy P.K.
 * Reference: Nature 275:673-674(1978).
 */
static const std::map<std::string, float> hphob_manavalan = {
	{"ALA", 12.97}, {"ARG", 11.72}, {"ASN", 11.42}, {"ASP", 10.85},
	{"CYS", 14.63}, {"GLN", 11.76}, {"GLU", 11.89}, {"GLY", 12.43},
	{"HIS", 12.16}, {"ILE", 15.67}, {"LEU", 14.90}, {"LYS", 11.36},
	{"MET", 14.39}, {"PHE", 14.00}, {"PRO", 11.37}, {"SER", 11.23},
	{"THR", 11.69}, {"TRP", 13.93}, {"TYR", 13.42}, {"VAL", 15.71}
};

/**
 * Hydrophobicity of physiological L-alpha amino acids
 * Author(s): Black S.D., Mould D.R.
 * Reference: Anal. Biochem. 193:72-82(1991).
 */
static const std::map<std::string, float> hphob_black = {
	{"ALA",  0.616}, {"ARG",  0.000}, {"ASN",  0.236}, {"ASP",  0.028},
	{"CYS",  0.680}, {"GLN",  0.251}, {"GLU",  0.043}, {"GLY",  0.501},
	{"HIS",  0.165}, {"ILE",  0.943}, {"LEU",  0.943}, {"LYS",  0.283},
	{"MET",  0.738}, {"PHE",  1.000}, {"PRO",  0.711}, {"SER",  0.359},
	{"THR",  0.450}, {"TRP",  0.878}, {"TYR",  0.880}, {"VAL",  0.825}
};

/**
 * Hydrophobicity scale (pi-r).
 * Author(s): Fauchere J.-L., Pliska V.E.
 * Reference: Eur. J. Med. Chem. 18:369-375(1983).
 */
static const std::map<std::string, float> hphob_fauchere = {
	{"ALA",  0.31}, {"ARG", -1.01}, {"ASN", -0.60}, {"ASP", -0.77},
	{"CYS",  1.54}, {"GLN", -0.22}, {"GLU", -0.64}, {"GLY",  0.00},
	{"HIS",  0.13}, {"ILE",  1.80}, {"LEU",  1.70}, {"LYS", -0.99},
	{"MET",  1.23}, {"PHE",  1.79}, {"PRO",  0.72}, {"SER", -0.04},
	{"THR",  0.26}, {"TRP",  2.25}, {"TYR",  0.96}, {"VAL",  1.22}
};

/**
 * Free energy of transfer from inside to outside of a globular protein.
 * This scale also provides information about the accessible and buried
 * amino acid residues of globular proteins.
 * Author(s): Janin J.
 * Reference: Nature 277:491-492(1979).
 */
static const std::map<std::string, float> hphob_janin = {
	{"ALA",  0.30}, {"ARG", -1.40}, {"ASN", -0.50}, {"ASP", -0.60},
	{"CYS",  0.90}, {"GLN", -0.70}, {"GLU", -0.70}, {"GLY",  0.30},
	{"HIS", -0.10}, {"ILE",  0.70}, {"LEU",  0.50}, {"LYS", -1.80},
	{"MET",  0.40}, {"PHE",  0.50}, {"PRO", -0.30}, {"SER", -0.10},
	{"THR", -0.20}, {"TRP",  0.30}, {"TYR", -0.40}, {"VAL",  0.60}
};

/**
 * Membrane buried helix parameter.
 * Author(s): Rao M.J.K., Argos P.
 * Reference: Biochim. Biophys. Acta 869:197-214(1986).
 */
static const std::map<std::string, float> hphob_rao = {
	{"ALA",  1.36}, {"ARG",  0.15}, {"ASN",  0.33}, {"ASP",  0.11},
	{"CYS",  1.27}, {"GLN",  0.33}, {"GLU",  0.25}, {"GLY",  1.09},
	{"HIS",  0.68}, {"ILE",  1.44}, {"LEU",  1.47}, {"LYS",  0.09},
	{"MET",  1.42}, {"PHE",  1.57}, {"PRO",  0.54}, {"SER",  0.97},
	{"THR",  1.08}, {"TRP",  1.00}, {"TYR",  0.83}, {"VAL",  1.37}
};

/**
 * Hydrophobicity scale (Contribution of hydrophobic interactions to the
 * stability of the globular conformation of proteins).
 * Author(s): Tanford C.
 * Reference: J. Am. Chem. Soc. 84:4240-4274(1962).
 */
static const std::map<std::string, float> hphob_tanford = {
	{"ALA",  0.62}, {"ARG", -2.53}, {"ASN", -0.78}, {"ASP", -0.09},
	{"CYS",  0.29}, {"GLN", -0.85}, {"GLU", -0.74}, {"GLY",  0.48},
	{"HIS", -0.40}, {"ILE",  1.38}, {"LEU",  1.53}, {"LYS", -1.50},
	{"MET",  0.64}, {"PHE",  1.19}, {"PRO",  0.12}, {"SER", -0.18},
	{"THR", -0.05}, {"TRP",  0.81}, {"TYR",  0.26}, {"VAL",  1.80}
};

/**
 * Antigenicity value X 10.
 * Welling et al. used information on the relative occurrence of amino acids
 * in antigenic regions to make a scale which is useful for prediction of
 * antigenic regions. This method is better than the Hopp-Woods scale of
 * hydrophobicity which is also used to identify antigenic regions.
 * Author(s): Welling G.W., Weijer W.J., Van der Zee R., Welling-Wester S.
 * Reference: FEBS Lett. 188:215-218(1985).
 */
static const std::map<std::string, float> hphob_welling = {
	{"ALA",  1.15}, {"ARG",  0.58}, {"ASN", -0.77}, {"ASP",  0.65},
	{"CYS", -1.20}, {"GLN", -0.11}, {"GLU", -0.71}, {"GLY", -1.84},
	{"HIS",  3.12}, {"ILE", -2.92}, {"LEU",  0.75}, {"LYS",  2.06},
	{"MET", -3.85}, {"PHE", -1.41}, {"PRO", -0.53}, {"SER", -0.26},
	{"THR", -0.45}, {"TRP", -1.14}, {"TYR",  0.13}, {"VAL", -0.13}
};

/**
 * Hydrophilicity scale derived from HPLC peptide retention times.
 * Author(s): Parker J.M.R., Guo D., Hodges R.S.
 * Reference: Biochemistry 25:5425-5431(1986).
 */
static const std::map<std::string, float> hphob_parker = {
	{"ALA",  2.10}, {"ARG",  4.20}, {"ASN",  7.00}, {"ASP", 10.00},
	{"CYS",  1.40}, {"GLN",  6.00}, {"GLU",  7.80}, {"GLY",  5.70},
	{"HIS",  2.10}, {"ILE", -8.00}, {"LEU", -9.20}, {"LYS",  5.70},
	{"MET", -4.20}, {"PHE", -9.20}, {"PRO",  2.10}, {"SER",  6.50},
	{"THR",  5.20}, {"TRP", -10.00}, {"TYR", -1.90}, {"VAL", -3.70}
};

/**
 * Hydrophobicity indices at ph 7.5 determined by HPLC.
 * Author(s): Cowan R., Whittaker R.G.
 * Reference: Peptide Research 3:75-80(1990).
 */
static const std::map<std::string, float> hphob_cowan2 = {
	{"ALA",  0.35}, {"ARG", -1.50}, {"ASN", -0.99}, {"ASP", -2.15},
	{"CYS",  0.76}, {"GLN", -0.93}, {"GLU", -1.95}, {"GLY",  0.00},
	{"HIS", -0.65}, {"ILE",  1.83}, {"LEU",  1.80}, {"LYS", -1.54},
	{"MET",  1.10}, {"PHE",  1.69}, {"PRO",  0.84}, {"SER", -0.63},
	{"THR", -0.27}, {"TRP",  1.35}, {"TYR",  0.39}, {"VAL",  1.32}
};
/**
 * The Engelman hydrophobicity scale, also known as the GES-scale,
 * is another scale which can be used for prediction of protein
 * hydrophobicity. As the Kyte-Doolittle scale, this scale is useful
 * for predicting transmembrane regions in proteins.
 * Author(s): Engelman, D. M., Steitz, T. A., and Goldman, A.
 * Reference: Annu Rev Biophys Biophys Chem, 15:321-353(1986).
 */
static const std::map<std::string, float> hphob_engelman = {
	{"ALA",	 1.60}, {"ARG",-12.30}, {"ASN", -4.80}, {"ASP", -9.20},
	{"CYS",	 2.00}, {"GLN",	-4.10}, {"GLU",	-8.20}, {"GLY",	 1.00},
	{"HIS",	-3.00}, {"ILE",	 3.10}, {"LEU",	 2.80}, {"LYS",	-8.80},
	{"MET",	 3.40}, {"PHE",	 3.70}, {"PRO",	-0.20}, {"SER",	 0.60},
	{"THR",	 1.20}, {"TRP",	 1.90}, {"TYR",	-0.70}, {"VAL",	 2.60}
};

/**
 * An optimal hydrophobicity scale based on 28 published scales.
 * Cornette et al. computed an optimal hydrophobicity scale based on 28
 * published scales. This optimized scale is also suitable for prediction
 * of alpha-helices in proteins.
 * Author(s): Cornette J.L., Cease K.B., Margalit H., Spouge J.L.,
 * Berzofsky J.A., DeLisi C.
 * Reference: J Mol Biol 1987;195:659-685.
 */
static const std::map<std::string, float> hphob_cornette = {
	{"ALA",  0.20},  {"ARG",  1.40},  {"ASN", -0.50},  {"ASP", -3.10},
	{"CYS",  4.10},  {"GLN", -2.80},  {"GLU", -1.80},  {"GLY",  0.00},
	{"HIS",  0.50},  {"ILE",  4.80},  {"LEU",  5.70},  {"LYS", -3.10},
	{"MET",  4.20},  {"PHE",  4.40},  {"PRO", -2.20},  {"SER", -0.50},
	{"THR", -1.90},  {"TRP",  1.00},  {"TYR",  3.20},  {"VAL",  4.70}
};

/**
 * Experimentally determined hydrophobicity scale for proteins
 * at membrane interfaces.
 * Author(s): Wimley W.C., White S.H.
 * Reference: Nature Struct Biol 1996;3:842-848.
 */
static const std::map<std::string, float> hphob_wimley = {
	{"ALA", -0.17}, {"ARG", -0.81}, {"ASN", -0.42}, {"ASP", -1.23},
	{"CYS",  0.24}, {"GLN", -0.58}, {"GLU", -2.02}, {"GLY", -0.01},
	{"HIS", -0.96}, {"ILE",  0.31}, {"LEU",  0.56}, {"LYS", -0.99},
	{"MET",  0.23}, {"PHE",  1.13}, {"PRO", -0.45}, {"SER", -0.13},
	{"THR", -0.14}, {"TRP",  1.85}, {"TYR",  0.94}, {"VAL", -0.07}
};

/**
 * Recognition of transmembrane helices by the endoplasmic reticulum translocon.
 * Author(s): Hessa T., Kim H., Bihlmaier K., Lundin C., Boekel J.,
 * Andersson H., Nilsson I., White SH., von Heijne G.
 * Reference: Nature 2005;433:377-381, supplementary data.
 */
static const std::map<std::string, float> hphob_hessa = {
	{"ALA",  0.11}, {"ARG",  2.58}, {"ASN",  2.05}, {"ASP",  3.49},
	{"CYS", -0.13}, {"GLN",  2.36}, {"GLU",  2.68}, {"GLY",  0.74},
	{"HIS",  2.06}, {"ILE", -0.60}, {"LEU", -0.55}, {"LYS",  2.71},
	{"MET", -0.10}, {"PHE", -0.32}, {"PRO",  2.23}, {"SER",  0.84},
	{"THR",  0.52}, {"TRP",  0.30}, {"TYR",  0.68}, {"VAL", -0.31}
};

/**
 * Binary atomic hydrophobicity values for each atom in the 20 natural amino acid residues.
 *
 * Each atom in a protein is identified as hydrophobic or hydrophilic. Atoms with a partial
 * charge magnitude ranging from 0 to 0.250 were considered non-polar and assigned a
 * hydrophobicity value of −1. Atoms with a partial charge magnitude greater than 0.250
 * were considered polar and were assigned a hydrophobicity value of +1.
 *
 * Hydrophobicity values are reported for each atom in the 20 natural amino acid residues.
 *
 * Lauren H. Kapcha, Peter J. Rossky, A Simple Atomic-Level Hydrophobicity Scale Reveals
 * Protein Interfacial Structure, Journal of Molecular Biology, Volume 426, Issue 2,
 * 23 January 2014, Pages 484-498, ISSN 0022-2836,
 * http://dx.doi.org/10.1016/j.jmb.2013.09.039.
 */

static const std::unordered_map<pair<string,string>, int> hphob_kapcha = {
		{ make_pair("ALA", "N"  ), +1 }, { make_pair("ALA", "CA" ), -1 }, { make_pair("ALA", "C"  ), +1 },
		{ make_pair("ALA", "O"  ), +1 }, { make_pair("ALA", "CB" ), -1 }, { make_pair("ALA", "H"  ), +1 },
		{ make_pair("ALA", "HA" ), -1 }, { make_pair("ALA", "HB1"), -1 }, { make_pair("ALA", "HB2"), -1 },
		{ make_pair("ALA", "HB3"), -1 },

		{ make_pair("ARG", "N"   ), +1 }, { make_pair("ARG", "CA"  ), -1 }, { make_pair("ARG", "C"   ), +1 },
		{ make_pair("ARG", "O"   ), +1 }, { make_pair("ARG", "CB"  ), -1 }, { make_pair("ARG", "CG"  ), -1 },
		{ make_pair("ARG", "CD"  ), -1 }, { make_pair("ARG", "NE"  ), +1 }, { make_pair("ARG", "CZ"  ), +1 },
		{ make_pair("ARG", "NH1" ), +1 }, { make_pair("ARG", "NH2" ), +1 }, { make_pair("ARG", "H"   ), +1 },
		{ make_pair("ARG", "HA"  ), -1 }, { make_pair("ARG", "HB2" ), -1 }, { make_pair("ARG", "HB3" ), -1 },
		{ make_pair("ARG", "HG2" ), -1 }, { make_pair("ARG", "HG3" ), -1 }, { make_pair("ARG", "HD2" ), -1 },
		{ make_pair("ARG", "HD3" ), -1 }, { make_pair("ARG", "HE"  ), +1 }, { make_pair("ARG", "HH11"), +1 },
		{ make_pair("ARG", "HH12"), +1 }, { make_pair("ARG", "HH21"), +1 }, { make_pair("ARG", "HH22"), +1 },

		{ make_pair("ASN", "N"   ), +1 }, { make_pair("ASN", "CA"  ), -1 }, { make_pair("ASN", "C"  ), +1 },
		{ make_pair("ASN", "O"   ), +1 }, { make_pair("ASN", "CB"  ), -1 }, { make_pair("ASN", "CG" ), +1 },
		{ make_pair("ASN", "OD1" ), +1 }, { make_pair("ASN", "ND2" ), +1 }, { make_pair("ASN", "H"  ), +1 },
		{ make_pair("ASN", "HA"  ), -1 }, { make_pair("ASN", "HB2" ), -1 }, { make_pair("ASN", "HB3"), -1 },
		{ make_pair("ASN", "HD21"), +1 }, { make_pair("ASN", "HD22"), +1 },

		{ make_pair("ASP", "N"  ), +1 }, { make_pair("ASP", "CA" ), -1 }, { make_pair("ASP", "C"  ), +1 },
		{ make_pair("ASP", "O"  ), +1 }, { make_pair("ASP", "CB" ), -1 }, { make_pair("ASP", "CG" ), +1 },
		{ make_pair("ASP", "OD1"), +1 }, { make_pair("ASP", "OD2"), +1 }, { make_pair("ASP", "H"  ), +1 },
		{ make_pair("ASP", "HA" ), -1 }, { make_pair("ASP", "HB2"), -1 }, { make_pair("ASP", "HB3"), -1 },

		{ make_pair("CYS", "N"  ), +1 }, { make_pair("CYS", "CA" ), -1 }, { make_pair("CYS", "C"  ), +1 },
		{ make_pair("CYS", "O"  ), +1 }, { make_pair("CYS", "CB" ), -1 }, { make_pair("CYS", "SG" ), +1 },
		{ make_pair("CYS", "H"  ), +1 }, { make_pair("CYS", "HA" ), -1 }, { make_pair("CYS", "HB2"), -1 },
		{ make_pair("CYS", "HB3"), -1 }, { make_pair("CYS", "HG" ), +1 },

		{ make_pair("GLN", "N"   ), +1 }, { make_pair("GLN", "CA"  ), -1 }, { make_pair("GLN", "C"   ), +1 },
		{ make_pair("GLN", "O"   ), +1 }, { make_pair("GLN", "CB"  ), -1 }, { make_pair("GLN", "CG"  ), -1 },
		{ make_pair("GLN", "CD"  ), +1 }, { make_pair("GLN", "OE1" ), +1 }, { make_pair("GLN", "NE2" ), +1 },
		{ make_pair("GLN", "H"   ), +1 }, { make_pair("GLN", "HA"  ), -1 }, { make_pair("GLN", "HB2" ), -1 },
		{ make_pair("GLN", "HB3" ), -1 }, { make_pair("GLN", "HG2" ), -1 }, { make_pair("GLN", "HG3" ), -1 },
		{ make_pair("GLN", "HE21"), +1 }, { make_pair("GLN", "HE22"), +1 },

		{ make_pair("GLU", "N"  ), +1 }, { make_pair("GLU", "CA" ), -1 }, { make_pair("GLU", "C"  ), +1 },
		{ make_pair("GLU", "O"  ), +1 }, { make_pair("GLU", "CB" ), -1 }, { make_pair("GLU", "CG" ), -1 },
		{ make_pair("GLU", "CD" ), +1 }, { make_pair("GLU", "OE1"), +1 }, { make_pair("GLU", "OE2"), +1 },
		{ make_pair("GLU", "H"  ), +1 }, { make_pair("GLU", "HA" ), -1 }, { make_pair("GLU", "HB2"), -1 },
		{ make_pair("GLU", "HB3"), -1 }, { make_pair("GLU", "HG2"), -1 }, { make_pair("GLU", "HG3"), -1 },

		{ make_pair("GLY", "N"  ), +1 }, { make_pair("GLY", "CA" ), -1 }, { make_pair("GLY", "C"  ), +1 },
		{ make_pair("GLY", "O"  ), +1 }, { make_pair("GLY", "H"  ), +1 }, { make_pair("GLY", "HA2"), -1 },
		{ make_pair("GLY", "HA3"), -1 },

		{ make_pair("HIS", "N"  ), +1 }, { make_pair("HIS", "CA" ), -1 }, { make_pair("HIS", "C"  ), +1 },
		{ make_pair("HIS", "O"  ), +1 }, { make_pair("HIS", "CB" ), -1 }, { make_pair("HIS", "CG" ), -1 },
		{ make_pair("HIS", "ND1"), +1 }, { make_pair("HIS", "CD2"), -1 }, { make_pair("HIS", "CE1"), +1 },
		{ make_pair("HIS", "NE2"), +1 }, { make_pair("HIS", "H"  ), +1 }, { make_pair("HIS", "HA" ), -1 },
		{ make_pair("HIS", "HB2"), -1 }, { make_pair("HIS", "HB3"), -1 }, { make_pair("HIS", "HD1"), +1 },
		{ make_pair("HIS", "HD2"), -1 }, { make_pair("HIS", "HE1"), +1 },

		{ make_pair("ILE", "CG1" ), -1 }, { make_pair("ILE", "CG2" ), -1 }, { make_pair("ILE", "CD1" ), -1 },
		{ make_pair("ILE", "H"   ), +1 }, { make_pair("ILE", "HA"  ), -1 }, { make_pair("ILE", "HB"  ), -1 },
		{ make_pair("ILE", "HG12"), -1 }, { make_pair("ILE", "HG13"), -1 }, { make_pair("ILE", "HG21"), -1 },
		{ make_pair("ILE", "HG22"), -1 }, { make_pair("ILE", "HG23"), -1 }, { make_pair("ILE", "HD11"), -1 },
		{ make_pair("ILE", "HD12"), -1 }, { make_pair("ILE", "HD13"), -1 },

		{ make_pair("LEU", "N"   ), +1 }, { make_pair("LEU", "CA"  ), -1 }, { make_pair("LEU", "C"   ), +1 },
		{ make_pair("LEU", "O"   ), +1 }, { make_pair("LEU", "CB"  ), -1 }, { make_pair("LEU", "CG"  ), -1 },
		{ make_pair("LEU", "CD1" ), -1 }, { make_pair("LEU", "CD2" ), -1 }, { make_pair("LEU", "H"   ), +1 },
		{ make_pair("LEU", "HA"  ), -1 }, { make_pair("LEU", "HB2" ), -1 }, { make_pair("LEU", "HB3" ), -1 },
		{ make_pair("LEU", "HG"  ), -1 }, { make_pair("LEU", "HD11"), -1 }, { make_pair("LEU", "HD12"), -1 },
		{ make_pair("LEU", "HD13"), -1 }, { make_pair("LEU", "HD21"), -1 }, { make_pair("LEU", "HD22"), -1 },
		{ make_pair("LEU", "HD23"), -1 },

		{ make_pair("LYS", "N"  ), +1 }, { make_pair("LYS", "CA" ), -1 }, { make_pair("LYS", "C"  ), +1 },
		{ make_pair("LYS", "O"  ), +1 }, { make_pair("LYS", "CB" ), -1 }, { make_pair("LYS", "CG" ), -1 },
		{ make_pair("LYS", "CD" ), -1 }, { make_pair("LYS", "CE" ), -1 }, { make_pair("LYS", "NZ" ), +1 },
		{ make_pair("LYS", "H"  ), +1 }, { make_pair("LYS", "HA" ), -1 }, { make_pair("LYS", "HB2"), -1 },
		{ make_pair("LYS", "HB3"), -1 }, { make_pair("LYS", "HG2"), -1 }, { make_pair("LYS", "HG3"), -1 },
		{ make_pair("LYS", "HD2"), -1 }, { make_pair("LYS", "HD3"), -1 }, { make_pair("LYS", "HE2"), -1 },
		{ make_pair("LYS", "HE3"), -1 }, { make_pair("LYS", "HZ1"), +1 }, { make_pair("LYS", "HZ2"), +1 },
		{ make_pair("LYS", "HZ3"), +1 },

		{ make_pair("MET", "N"  ), +1 }, { make_pair("MET", "CA" ), -1 }, { make_pair("MET", "C"  ), +1 },
		{ make_pair("MET", "O"  ), +1 }, { make_pair("MET", "CB" ), -1 }, { make_pair("MET", "CG" ), -1 },
		{ make_pair("MET", "SD" ), +1 }, { make_pair("MET", "CE" ), -1 }, { make_pair("MET", "H"  ), +1 },
		{ make_pair("MET", "HA" ), -1 }, { make_pair("MET", "HB2"), -1 }, { make_pair("MET", "HB3"), -1 },
		{ make_pair("MET", "HG2"), -1 }, { make_pair("MET", "HG3"), -1 }, { make_pair("MET", "HE1"), -1 },
		{ make_pair("MET", "HE2"), -1 }, { make_pair("MET", "HE3"), -1 },

		{ make_pair("PHE", "N"  ), +1 }, { make_pair("PHE", "CA" ), -1 }, { make_pair("PHE", "C"  ), +1 },
		{ make_pair("PHE", "O"  ), +1 }, { make_pair("PHE", "CB" ), -1 }, { make_pair("PHE", "CG" ), -1 },
		{ make_pair("PHE", "CD1"), -1 }, { make_pair("PHE", "CD2"), -1 }, { make_pair("PHE", "CE1"), -1 },
		{ make_pair("PHE", "CE2"), -1 }, { make_pair("PHE", "CZ" ), -1 }, { make_pair("PHE", "H"  ), +1 },
		{ make_pair("PHE", "HA" ), -1 }, { make_pair("PHE", "HB2"), -1 }, { make_pair("PHE", "HB3"), -1 },
		{ make_pair("PHE", "HD1"), -1 }, { make_pair("PHE", "HD2"), -1 }, { make_pair("PHE", "HE1"), -1 },
		{ make_pair("PHE", "HE2"), -1 }, { make_pair("PHE", "HZ" ), -1 },

		{ make_pair("PRO", "N"  ), -1 }, { make_pair("PRO", "CA" ), -1 }, { make_pair("PRO", "C"  ), +1 },
		{ make_pair("PRO", "O"  ), +1 }, { make_pair("PRO", "CB" ), -1 }, { make_pair("PRO", "CG" ), -1 },
		{ make_pair("PRO", "CD" ), -1 }, { make_pair("PRO", "HA" ), -1 }, { make_pair("PRO", "HB2"), -1 },
		{ make_pair("PRO", "HB3"), -1 }, { make_pair("PRO", "HG2"), -1 }, { make_pair("PRO", "HG3"), -1 },
		{ make_pair("PRO", "HD2"), -1 }, { make_pair("PRO", "HD3"), -1 },

		{ make_pair("SER", "N"  ), +1 }, { make_pair("SER", "CA" ), -1 }, { make_pair("SER", "C"  ), +1 },
		{ make_pair("SER", "O"  ), +1 }, { make_pair("SER", "CB" ), -1 }, { make_pair("SER", "OG" ), +1 },
		{ make_pair("SER", "H"  ), +1 }, { make_pair("SER", "HA" ), -1 }, { make_pair("SER", "HB2"), -1 },
		{ make_pair("SER", "HB3"), -1 }, { make_pair("SER", "HG" ), +1 },

		{ make_pair("THR", "N"   ), +1 }, { make_pair("THR", "CA"  ), -1 }, { make_pair("THR", "C"   ), +1 },
		{ make_pair("THR", "O"   ), +1 }, { make_pair("THR", "CB"  ), -1 }, { make_pair("THR", "OG1" ), +1 },
		{ make_pair("THR", "CG2" ), -1 }, { make_pair("THR", "H"   ), +1 }, { make_pair("THR", "HA"  ), -1 },
		{ make_pair("THR", "HB"  ), -1 }, { make_pair("THR", "HG1" ), +1 }, { make_pair("THR", "HG21"), -1 },
		{ make_pair("THR", "HG22"), -1 }, { make_pair("THR", "HG23"), -1 },

		{ make_pair("TRP", "N"  ), +1 }, { make_pair("TRP", "CA" ), -1 }, { make_pair("TRP", "C"  ), +1 },
		{ make_pair("TRP", "O"  ), +1 }, { make_pair("TRP", "CB" ), -1 }, { make_pair("TRP", "CG" ), -1 },
		{ make_pair("TRP", "CD1"), -1 }, { make_pair("TRP", "CD2"), -1 }, { make_pair("TRP", "NE1"), +1 },
		{ make_pair("TRP", "CE2"), -1 }, { make_pair("TRP", "CE3"), -1 }, { make_pair("TRP", "CZ2"), -1 },
		{ make_pair("TRP", "CZ3"), -1 }, { make_pair("TRP", "CH2"), -1 }, { make_pair("TRP", "H"  ), +1 },
		{ make_pair("TRP", "HA" ), -1 }, { make_pair("TRP", "HB2"), -1 }, { make_pair("TRP", "HB3"), -1 },
		{ make_pair("TRP", "HD1"), -1 }, { make_pair("TRP", "HE1"), +1 }, { make_pair("TRP", "HE3"), -1 },
		{ make_pair("TRP", "HZ2"), -1 }, { make_pair("TRP", "HZ3"), -1 },

		{ make_pair("TYR", "N"  ), +1 }, { make_pair("TYR", "CA" ), -1 }, { make_pair("TYR", "C"  ), +1 },
		{ make_pair("TYR", "O"  ), +1 }, { make_pair("TYR", "CB" ), -1 }, { make_pair("TYR", "CG" ), -1 },
		{ make_pair("TYR", "CD1"), -1 }, { make_pair("TYR", "CD2"), -1 }, { make_pair("TYR", "CE1"), -1 },
		{ make_pair("TYR", "CE2"), -1 }, { make_pair("TYR", "CZ" ), -1 }, { make_pair("TYR", "OH" ), +1 },
		{ make_pair("TYR", "H"  ), +1 }, { make_pair("TYR", "HA" ), -1 }, { make_pair("TYR", "HB2"), -1 },
		{ make_pair("TYR", "HB3"), -1 }, { make_pair("TYR", "HD1"), -1 }, { make_pair("TYR", "HD2"), -1 },
		{ make_pair("TYR", "HE1"), -1 }, { make_pair("TYR", "HE2"), -1 }, { make_pair("TYR", "HH" ), +1 },

		{ make_pair("VAL", "N"   ), +1 }, { make_pair("VAL", "CA"  ), -1 }, { make_pair("VAL", "C"   ), +1 },
		{ make_pair("VAL", "O"   ), +1 }, { make_pair("VAL", "CB"  ), -1 }, { make_pair("VAL", "CG1" ), -1 },
		{ make_pair("VAL", "CG2" ), -1 }, { make_pair("VAL", "H"   ), +1 }, { make_pair("VAL", "HA"  ), -1 },
		{ make_pair("VAL", "HB"  ), -1 }, { make_pair("VAL", "HG11"), -1 }, { make_pair("VAL", "HG12"), -1 },
		{ make_pair("VAL", "HG13"), -1 }, { make_pair("VAL", "HG21"), -1 }, { make_pair("VAL", "HG22"), -1 },
		{ make_pair("VAL", "HG23"), -1 }
};
#endif /* HYDROPHOBICITY_H_ */
