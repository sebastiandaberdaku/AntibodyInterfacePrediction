/*
 * weights.h
 *
 *  Created on: Jul 5, 2016
 *      Author: sebastian
 */

#ifndef WEIGHTS_H_
#define WEIGHTS_H_

#include <vector>
#include <fstream>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/filesystem.hpp>
#include <string>
#include <iostream>
#include <sstream>

#include "DockingMethods/SurfacePair.h"


#include "utils/progressBar.h"


using namespace std;


void export_descriptors(vector<CompactPatchDescriptor> & descriptors, size_t order, string const & outname) {
	ofstream out_descriptors;
	out_descriptors.open("./" + outname + "_train_descriptors_N" + to_string(order) + ".txt", ofstream::out | ofstream::trunc);
	for (auto & d : descriptors) {
		CompactPatchDescriptor::to_ostream(out_descriptors,  d, order);
			out_descriptors << endl;
	}
	out_descriptors.close();

}




#endif /* WEIGHTS_H_ */
