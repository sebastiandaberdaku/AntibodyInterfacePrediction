/*
 * hash.h
 *
 *  Created on: Mar 3, 2016
 *      Author: sebastian
 */

#ifndef HASH_H_
#define HASH_H_



namespace std {
template<>
struct hash<pair<string, string>> {
	size_t operator()(pair<string,string> const & p) const {
		size_t result = hash<string>()(p.first);
		boost::hash_combine(result, p.second);
		return result;
	}
};

template<>
struct equal_to<pair<string, string>> {
	bool operator()(pair<string, string> const & p, pair<string, string> const & q) const {
		return (p.first == q.first) && (p.second == q.second);
	}
};
}

#endif /* HASH_H_ */
