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


void export_descriptors(vector<CompactPatchDescriptor> & descriptors, size_t order, string const & outname) {
	ofstream out_descriptors, out_centers, out_truth;
	out_descriptors.open("./" + outname + "_test_descriptors_N" + to_string(order) + ".txt", ofstream::out | ofstream::trunc);
	out_centers.open("./" + outname + "_patch_centers.txt", ofstream::out | ofstream::trunc);
	out_truth.open("./" + outname + "_patch_truth.txt", ofstream::out | ofstream::trunc);

	for (auto & d : descriptors) {
		CompactPatchDescriptor::to_ostream(out_descriptors,  d, order);
			out_descriptors << endl;
		out_centers << d.center.x << "\t" << d.center.y << "\t" << d.center.z << "\t" << endl;
		if (d.isInterface)
			out_truth << "+1" << endl;
		else
			out_truth << "-1" << endl;
	}
	out_descriptors.close();
	out_centers.close();
	out_truth.close();

}




#endif /* WEIGHTS_H_ */
