/*
 * CustomValidators.h
 *
 *  Created on: Jan 11, 2014
 *      Author: sebastian
 */
/**
 * Here we overload the validate function in order to accept only valid
 * command line input parameters.
 *
 * For a detailed description on overloading the default validator see:
 * http://www.boost.org/doc/libs/1_55_0/doc/html/program_options/howto.html#idp163429032
 */
#ifndef CUSTOMVALIDATORS_H_
#define CUSTOMVALIDATORS_H_

#include <boost/program_options.hpp>
#include <boost/regex.hpp>
#include <string>
using namespace std;
using namespace boost;
using namespace boost::program_options;

struct input_PDB_PQR_filename {
	input_PDB_PQR_filename(string const & val) :
			filename(val) {	}
	string filename;
	friend ostream& operator <<(ostream& s, const input_PDB_PQR_filename& ipqrf) {
		s << ipqrf.filename;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * PQR input filenames.
 *
 * The function takes four parameters. The first is the storage for
 * the value, and in this case is either empty or contains an instance
 * of the magic_number class. The second is the list of strings found
 * in the next occurrence of the option. The remaining two parameters
 * are needed to workaround the lack of partial template specialization
 * and partial function template ordering on some compilers.
 *
 * The function first checks that we don't try to assign to the same
 * option twice. Then it checks that only a single string was passed
 * in. Next the string is verified. If that test is passed, the parsed
 * value is stored into the v variable.
 */

void validate(any& v, vector<string> const& values,
		input_PDB_PQR_filename* /* target_type */, int) {
	static regex r("[^[.NUL.]]+(?:\\.pqr|\\.PQR|\\.pdb|\\.PDB)");
	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);
	// Do regex match
	if (regex_match(s, r))
		v = any(input_PDB_PQR_filename(s));
	 else
		throw validation_error(validation_error::invalid_option_value);
};
struct filename {
	filename(string const & outname) :
		fname(outname) {	}
	string fname;
	friend ostream& operator <<(ostream& s, const filename& out) {
		s << out.fname;
		return s;
	}
};
void validate(any& v, vector<string> const& values,
		filename* /* target_type */, int) {
	static regex r("[^[.NUL.]]+");
	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);
	// Do regex match
	if (regex_match(s, r))
		v = any(filename(s));
	 else
		throw validation_error(validation_error::invalid_option_value);
};

struct input_DX_filename {
	input_DX_filename(string const& val) :
			filename(val) {	}
	string filename;
	friend ostream& operator <<(ostream& s, const input_DX_filename& idxf) {
		s << idxf.filename;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * openDX input filenames.
 *
 * The function takes four parameters. The first is the storage for
 * the value, and in this case is either empty or contains an instance
 * of the magic_number class. The second is the list of strings found
 * in the next occurrence of the option. The remaining two parameters
 * are needed to workaround the lack of partial template specialization
 * and partial function template ordering on some compilers.
 *
 * The function first checks that we don't try to assign to the same
 * option twice. Then it checks that only a single string was passed
 * in. Next the string is verified. If that test is passed, the parsed
 * value is stored into the v variable.
 */
void validate(any& v, vector<string> const& values,
		input_DX_filename* /* target_type */, int) {
	static regex r("[^[.NUL.]]+\\.dx");
	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);
	// Do regex match
	if (regex_match(s, r))
		v = any(input_DX_filename(s));
	 else
		throw validation_error(validation_error::invalid_option_value);
};

struct max_order {
public:
	max_order(int n) :
			n(n) { }
	int n;
	friend ostream& operator <<(ostream& s, const max_order& mo) {
		s << mo.n;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * maximum Zernike descriptor order input parameters.
 */
void validate(any& v, vector<string> const& values,
		max_order* /* target_type */, int) {
	static regex r("[0-9]+");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		int temp = lexical_cast<int>(s);
		if (temp > 0)
			v = any(max_order(lexical_cast<int>(s)));
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

struct dockingPoses_param {
public:
	dockingPoses_param(int n) :
			n(n) { }
	int n;
	friend ostream& operator <<(ostream& s, const dockingPoses_param& mo) {
		s << mo.n;
		return s;
	}
};
//TODO

/**
 * Here we overload the validate function in order to accept only valid
 * number of docking poses input parameters.
 */
void validate(any& v, vector<string> const& values,
		dockingPoses_param* /* target_type */, int) {
	static regex r("[0-9]+");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		int temp = lexical_cast<int>(s);
		if (temp > 0)
			v = any(dockingPoses_param(lexical_cast<int>(s)));
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

struct bestPairs_param {
public:
	bestPairs_param(float p) :
			p(p) { }
	float p;
	friend ostream& operator <<(ostream& s, const bestPairs_param& bp) {
		s << bp.p;
		return s;
	}
};
//TODO
/**
 * Here we overload the validate function in order to accept only valid
 * best-matching surface patch descriptors input parameters.
 */
void validate(any& v, vector<string> const& values,
		bestPairs_param* /* target_type */, int) {
	static regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		int temp = lexical_cast<float>(s);
		if (temp > 0 && temp <= 1.0)
			v = any(bestPairs_param(lexical_cast<int>(s)));
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

struct probe_radius {
public:
	probe_radius(float p) :
			p(p) { }
	float p;
	friend ostream& operator <<(ostream& s, const probe_radius& pr) {
		s << pr.p;
		return s;
	}
};
/**
 * Here we overload the validate function in order to accept only valid
 * probe radius parameters.
 */

void validate(any& v, vector<string> const& values,
		probe_radius* /* target_type */, int) {
	static regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		float temp = lexical_cast<float>(s);
		if (temp > 0 && temp <= 2.0)
			v = any(probe_radius(temp));
		else
			throw validation_error(validation_error::invalid_option_value);
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

struct patch_rd {
public:
	patch_rd(float p) :
			p(p) { }
	float p;
	friend ostream& operator <<(ostream& s, const patch_rd& pr) {
		s << pr.p;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * patch radius and minimum distance between patch centers parameters.
 */

void validate(any& v, vector<string> const& values,
		patch_rd* /* target_type */, int) {
	static regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		float temp = lexical_cast<float>(s);
		if (temp > 0)
			v = any(patch_rd(temp));
		else
			throw validation_error(validation_error::invalid_option_value);
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

struct patch_threshold {
public:
	patch_threshold(float t) :
			t(t) { }
	float t;
	friend ostream& operator <<(ostream& s, const patch_threshold& pt) {
		s << pt.t;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * patch radius and minimum distance between patch centers parameters.
 */

void validate(any& v, vector<string> const& values,
		patch_threshold* /* target_type */, int) {
	static regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		float temp = lexical_cast<float>(s);
		if (temp > 0 && temp <= 1)
			v = any(patch_threshold(temp));
		else
			throw validation_error(validation_error::invalid_option_value);
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};
struct resolution_param {
public:
	resolution_param(float r) :
			r(r) { }
	float r;
	friend ostream& operator <<(ostream& s, const resolution_param& rp) {
		s << rp.r;
		return s;
	}
};

/**
 * Here we overload the validate function in order to accept only valid
 * resolution values.
 */

void validate(any& v, vector<string> const& values,
		resolution_param* /* target_type */, int) {
	static regex r("[+-]?[0-9]*\\.?[0-9]*");

	// Make sure no previous assignment to 'v' was made.
	validators::check_first_occurrence(v);

	// Extract the first string from 'values'. If there is more than
	// one string, it's an error, and exception will be thrown.
	string const& s = validators::get_single_string(values);

	if (regex_match(s, r)) {
		float temp = lexical_cast<float>(s);
		if (temp > 0)
			v = any(resolution_param(temp));
		else
			throw validation_error(validation_error::invalid_option_value);
	}
	else
		throw validation_error(validation_error::invalid_option_value);
};

#endif /* CUSTOMVALIDATORS_H_ */
