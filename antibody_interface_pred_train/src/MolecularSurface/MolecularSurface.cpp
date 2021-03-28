/*
 * MolecularSurface.cpp
 *
 *  Created on: Dec 23, 2013
 *      Author: daberdaku
 */
/**
 * Implementation of the MolecularSurface class. This class
 * defines the molecular surface object, with all the methods
 * needed for its calculation.
 */
#include "../ASP/AtomicSolvationParameters.h"
#include "../exceptions/ParsingPQRException.h"
#include "../Geometry/HierarchicalQueue.h"
#include "../Geometry/normalizeGrid.h"
#include "../Geometry/partialEDMap.h"
#include "../Geometry/SeedFill3D/rapid3DCPKModel.h"
#include "../Geometry/SeedFill3D/rapid3DSurfaceExtract.h"
#include "../Geometry/Sphere/DrawSphere.h"
#include "../Geometry/Sphere/SphereOffsets.h"
#include "../MolecularSurface/MolecularSurface.h"
#include "../utils/makeDirectory.h"
#include "../utils/numerical_limits.h"
#include "../Zernike/ZernikeDescriptor.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>


#include "../Atom/AAindex1.h"

#ifdef PROGRESS_BAR
#include "../utils/progressBar.h"
#endif

/**
 * Constructor of the class. Initializes the voxel grid and other data structures.
 *
 * \param atoms				List containing the atoms of the molecule.
 * \param probeRadius		Desired probe radius of the solvent rolling sphere.
 * \param resolution		resolution^3 = number of voxels per Å^3
 * \param isPQR				true if the atom list comes from a PQR file, false
 * 							if it comes from a PDB file
 * \param addProbeRadius	true: add the probe radius to the atomic radius
 * 							values, SAS and MS (SES) have this value set to true;
 * 							false: use original atomic radius values.
 * \param length			Length of the voxel grid.
 * \param width				Width of the voxel grid.
 * \param height			Height of the voxel grid.
 * \param translation		Translation vector to the grid's coordinate system
 * 							(already scaled).
 *
 */
MolecularSurface::MolecularSurface(std::vector<atom> const & atoms,
		float probeRadius, float resolution,
		uint16_t length, uint16_t width, uint16_t height, point3D translation) :
		atoms(&atoms), probeRadius(probeRadius), resolution(resolution),
		plength(length), pwidth(width), pheight(height), ptran(translation),
		cpkModel(NULL), surface(NULL), DXPotentials(NULL) {
	if (atoms.empty())
		throw invalid_argument("MolecularSurface::MolecularSurface() - The molecule has no atoms!");
} /* MolecularSurface() */

/**
 * Destructor of the class.
 */
MolecularSurface::~MolecularSurface() {
	if (cpkModel != NULL)
		delete cpkModel;
	if (surface != NULL)
		delete surface;
	if (DXPotentials != NULL)
		delete DXPotentials;
} /* ~MolecularSurface() */

/** Calculates a bounding box containing the molecule. This method
 * calculates the maximal and minimal coordinates reached by the
 * molecule by checking all the coordinates of the atoms composing it.
 *
 * \param atomsInModel		List containing all the model's atoms.
 * \param probeRadius		Desired probe radius of the solvent rolling sphere.
 * \param resolution		resolution^3 = number of voxels per Å^3
 * \param length			(return value) the length of the voxel grid.
 * \param width				(return value) the width of the voxel grid.
 * \param height			(return value) the height of the voxel grid.
 * \param translation		(return value) the translation vector to the
 * 							grid's coordinate system (already scaled).
 * \param max_atm_radii		(return value) the maximum atomic radii in
 * 							atomsInModel
 */
void MolecularSurface::boundingBox(std::vector<atom> const & atomsInModel,
		float probeRadius, float patchRadius, float resolution, uint16_t & length,
		uint16_t & width, uint16_t & height, point3D & translation,
		float & max_atm_radius, float & min_atm_radius) {
	if (atomsInModel.empty())
		throw ParsingPQRException("No atoms in PQR model.",
				"MolecularSurface::boundingBox", "Incorrect PQR file format.");
	point3D minp(max_float, max_float, max_float);
	point3D maxp(min_float, min_float, min_float);
	for (std::vector<atom>::const_iterator atm = atomsInModel.begin();
			atm != atomsInModel.end(); ++atm) {
			if (atm->x < minp.x)
				minp.x = atm->x;
			if (atm->y < minp.y)
				minp.y = atm->y;
			if (atm->z < minp.z)
				minp.z = atm->z;
			if (atm->x > maxp.x)
				maxp.x = atm->x;
			if (atm->y > maxp.y)
				maxp.y = atm->y;
			if (atm->z > maxp.z)
				maxp.z = atm->z;
			if (max_atm_radius < atm->radius)
				max_atm_radius = atm->radius;
			if (min_atm_radius > atm->radius && atm->radius > 0)
				min_atm_radius = atm->radius;
	}
	float d = max_atm_radius + 2 * probeRadius + 0.5;
	maxp += point3D(d, d, d);
	minp -= point3D(d, d, d);

	/* transformation values */
	translation = -minp;

	/* bounding box dimensions */
	double boxLength = ceil(resolution * (maxp.x - minp.x));
	double boxWidth  = ceil(resolution * (maxp.y - minp.y));
	double boxHeight = ceil(resolution * (maxp.z - minp.z));

	if (boxLength <= UINT16_MAX && boxWidth <= UINT16_MAX && boxHeight <= UINT16_MAX) {
		length = static_cast<uint16_t>(boxLength);
		width  = static_cast<uint16_t>(boxWidth);
		height = static_cast<uint16_t>(boxHeight);
	} else {
		std::stringstream ss;
		ss << "MolecularSurface::boundingBox() - ";
		ss << "The bounding box's dimensions exceed the maximum value of " << UINT16_MAX << ". ";
		ss << "Try setting a lower \"resolution\" value.";
		throw std::invalid_argument(ss.str());
	}
} /* boundingBox() */

/** Fills the voxels in the grid occupied by the molecule (protein).
 * This method implements a space-filling algorithm which is the
 * preliminary step for our grid-based macro-molecular surface
 * generation.
 *
 * In chemistry, a space-filling model, also known as a calotte model,
 * is a type of three-dimensional molecular model where the atoms are
 * represented by spheres whose radii are proportional to the radii of
 * the atoms and whose center-to-center distances are proportional to
 * the distances between the atomic nuclei, all in the same scale.
 * Atoms of different chemical elements are usually represented by
 * spheres of different colors.
 *
 * Calotte models are distinguished from other 3D representations,
 * such as the ball-and-stick and skeletal models, by the use of "full
 * size" balls for the atoms. They are useful for visualizing the
 * effective shape and relative dimensions of the molecule, in particular
 * the region of space occupied by it. On the other hand, calotte models
 * do not show explicitly the chemical bonds between the atoms, nor the
 * structure of the molecule beyond the first layer of atoms.
 *
 * Space-filling models are also called CPK models after the chemists
 * Robert Corey, Linus Pauling and Walter Koltun, who pioneered their use.
 */
void MolecularSurface::createCPKModel() {
	if (cpkModel != NULL)
		delete cpkModel;
	cpkModel = new voxelGrid(plength, pwidth, pheight);
	uint16_t cx, cy, cz; /**< Discretized coordinates of the atom's center. */
	uint16_t radius; /**< radius */
	std::vector<atom>::const_iterator atm;
	/* For every atom in our list, calculate the voxels it occupies. */

	for (atm = atoms->begin(); atm != atoms->end(); ++atm) {
		/* Hydrogens sometimes have null radii in PQR files */
		if (atm->hasZeroRadius())
			continue;

		/* Translate and discretize the coordinates */
		cx = static_cast<uint16_t>((atm->x + ptran.x) * resolution + 0.5);
		cy = static_cast<uint16_t>((atm->y + ptran.y) * resolution + 0.5);
		cz = static_cast<uint16_t>((atm->z + ptran.z) * resolution + 0.5);

		radius = static_cast<uint16_t>((atm->radius + probeRadius) * resolution + 0.5);

		DrawSphere(*cpkModel, cx, cy, cz, radius);
	}
} /* createCPKModel() */
void MolecularSurface::createCPKModel_old() {
	if (cpkModel != NULL)
		delete cpkModel;

	cpkModel = new voxelGrid(plength, pwidth, pheight);
	uint16_t cx, cy, cz; /**< Discretized coordinates of the atom's center. */
	uint16_t radius; /**< radius */
	uint16_t y2_max; /**< maximum value for the squared y coordinate */
	uint16_t z2_max; /**< maximum value for the squared z coordinate */
	std::vector<atom>::const_iterator atm;
	/* For every atom in our list, calculate the voxels it occupies. */

	for (atm = atoms->begin(); atm != atoms->end(); ++atm) {
		/* Hydrogens sometimes have null radii in PQR files */
		if (atm->hasZeroRadius())
			continue;

		/* Translate and discretize the coordinates */
		cx = static_cast<uint16_t>((atm->x + ptran.x) * resolution + 0.5);
		cy = static_cast<uint16_t>((atm->y + ptran.y) * resolution + 0.5);
		cz = static_cast<uint16_t>((atm->z + ptran.z) * resolution + 0.5);

		radius = static_cast<uint16_t>((atm->radius + probeRadius) * resolution + 0.5);
		// x^2 + y^2 + z^2 = r^2
		for (uint16_t x = 0; x <= radius; ++x) {
			// y^2 + z^2 = r^2 - x^2
			y2_max = radius * radius - x * x;
			for (uint16_t y = 0; y * y <= y2_max; ++y) {
				// z^2 = r^2 - x^2 - y^2
				z2_max = radius * radius - x * x - y * y;
				for (uint16_t z = 0; z * z <= z2_max; ++z) {
					cpkModel->setVoxel(cx + x, cy + y, cz + z);
					cpkModel->setVoxel(cx + x, cy + y, cz - z);
					cpkModel->setVoxel(cx + x, cy - y, cz + z);
					cpkModel->setVoxel(cx + x, cy - y, cz - z);
					cpkModel->setVoxel(cx - x, cy + y, cz + z);
					cpkModel->setVoxel(cx - x, cy + y, cz - z);
					cpkModel->setVoxel(cx - x, cy - y, cz + z);
					cpkModel->setVoxel(cx - x, cy - y, cz - z);
				}
			}
		}
	}
} /* createCPKModel() */

/** Method that removes any spurious internal voids. Macromolecules
 * can have solvent-excluded cavities and voids, which might generate
 * spurious surfaces inside the real molecular surface. This method
 * calls a simple three-dimensional flood fill algorithm, which
 * "colors" the voxels on the outside of the most external protein
 * surface. The "coloring" process starts from one of the eight
 * vertices of the bounding box. Note that, if an adequately high
 * resolution is used, the voxels at the vertices of the bounding box
 * will never belong to any atom. Just think of fitting some spheres
 * inside a box: there will always be some free space at its edges.
 */
void MolecularSurface::fillInternalCavities_old() {
	assert(!cpkModel->getVoxel(0, 0, 0));

	floodFill3D(*cpkModel, voxel(0, 0, 0));
}
void MolecularSurface::fillInternalCavities() {
	/* Determine the starting point for the flood fill algorithm.
	 * Start by testing an edge voxel of the bounding box.
	 * If an adequate resolution is used, the edges should not be
	 * occupied by any atom. This possibility indicates an erroneous
	 * calculation of the protein surface. */
	assert(!cpkModel->getVoxel(0, 0, 0));

	voxelGrid temp(*cpkModel);
	/* Call to the three-dimensional flood fill algorithm to
	 * "color" the all voxels lying between the bounding box's
	 * boundaries and the outer surface of the volumetric model.
	 */
	rapid3DCPKModel::cpkModelFill(temp, *cpkModel, voxel(0, 0, 0));

} /* fillInternalCavities() */

/** Simple three-dimensional flood fill algorithm. The "coloring" starts
 * from startingVoxel, and proceeds by setting the occupied variable to
 * "false" for all the voxels lying outside the molecular surface. The
 * algorithm uses an auxiliary queue data structure for storing the next
 * voxels to be "colored".
 *
 * \param startingVoxel		The starting voxel for the algorithm.
 * \param grid		pointer to the voxelGrid
 * \param lenght	length of the voxel grid
 * \param width		width of the voxel grid
 * \param height	height of the voxel grid
 */
void MolecularSurface::floodFill3D(voxelGrid & grid, voxel const & startingVoxel) {
	/* Sets temp equal to grid. Both variables will describe
	 * the same volume. */
	voxelGrid temp(grid);

	/* initialize */
	grid.set_all();

	std::queue<voxel> voxelQueue;
	voxel currentVoxel = startingVoxel;

	voxel nb;
	voxelQueue.push(currentVoxel);
	grid.clearVoxel(currentVoxel.ix, currentVoxel.iy, currentVoxel.iz);

	while (!voxelQueue.empty()) {
		currentVoxel = voxelQueue.front();
		voxelQueue.pop();
		/* for all currentVoxel's neighbors */
		for (int i = 0; i < 6; ++i) {
			/* try to fill currentVoxel's neighbors */
			nb = currentVoxel + voxel_offset::neighbours[i];
			/* If not going out of bounds (i.e. if the neighbor exists),
			 * if the current neighbor has not been colored yet, and if
			 * it is not occupied by an atom, put it in the queue. */
			if (nb.iz < grid.height && nb.iy < grid.width	&& nb.ix < grid.length) {
				if (grid.getVoxel(nb) && !temp.getVoxel(nb)) {
					// color it!
					grid.clearVoxel(nb);
					voxelQueue.push(nb);
				} // if
			} // if
		} // for
	} // while
} /* floodFill3D() */

/**
 * Build the boundary of the solid created with the space-filling algorithm.
 * A voxel belongs to the boundary if it has at least one neighbor which is not
 * occupied by any atom.
 */
void MolecularSurface::buildSurface() {
	if (surface != NULL)
		delete surface;
	surface = new voxelGrid(plength, pwidth, pheight); /**< already initialized to all zeros */
	uint16_t cx, cy, cz;
	/*
	 * find a first voxel belonging to the CPKmodel
	 * and use it as a seed;
	 */
	voxelGrid temp(*cpkModel);
	/*
	 * the model could have disjoint parts,
	 * all atom centers must be used as seeds
	 */
	for (auto const & atm : *atoms) {
		/* Hydrogens usually have null radii in PQR files */
		if (atm.hasZeroRadius())
			continue;
		/* Translate and discretize the coordinates */
		cx = static_cast<uint16_t>((atm.x + ptran.x) * resolution + 0.5);
		cy = static_cast<uint16_t>((atm.y + ptran.y) * resolution + 0.5);
		cz = static_cast<uint16_t>((atm.z + ptran.z) * resolution + 0.5);
		if (temp.getVoxel(cx, cy, cz))
			rapid3DSurfaceExtract::surfaceExtract3D(*cpkModel, temp, *surface, voxel(cx, cy, cz));
	}

} /* buildSurface() */
void MolecularSurface::buildSurface_old() {
	if (surface != NULL)
		delete surface;
	surface = new voxelGrid(plength, pwidth, pheight);

	int ii; /**< Neighbor index - each voxel has at most 26 neighbors. */
	bool flagbound; /**< Boundary switch - if an internal voxel has an external
	 	 	 	 	   * neighbor than stop checking the other neighbors. */
	voxel nb;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				/* if the current voxel is inside the molecule */
				if (cpkModel->getVoxel(i, j, k)) {
					flagbound = false;
					ii = 0;
					/* while not a boundary voxel, and while not all of its
					 * neighbors have been checked */
					while (!flagbound && ii < 26) {
						nb = voxel(i, j, k) + voxel_offset::neighbours[ii];
						if (!cpkModel->getVoxel(nb)) {
							surface->setVoxel(i, j, k); /* is a boundary voxel */
							/* stop checking the remaining neighbors */
							flagbound = true;
						} // if
						else {
							++ii; /* go on with the next neighbor */
						} // else
					} // while
				} // if
			} // for k
		} // for j
	} // for i
} /* buildSurface() */
/**
 * Very similar to the regionGrowingEDT() method, but instead
 * of employing the partialEDMap data structure, it uses a 3D
 * array of voxels to keep track of the nearest boundary point
 * of each voxel, and a 3D array of integers to keep track of
 * the squared distance from the nearest boundary point for
 * each voxel.
 */
void MolecularSurface::fastRegionGrowingEDT() {
	voxel ***nbv = new voxel**[plength];
	for (int i = 0; i < plength; ++i) {
		nbv[i] = new voxel*[pwidth];
	}
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			nbv[i][j] = new voxel[pheight];
		}
	}
	uint16_t ***dMap = new uint16_t**[plength];
	for (int i = 0; i < plength; ++i) {
		dMap[i] = new uint16_t*[pwidth];
	}
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			dMap[i][j] = new uint16_t[pheight];
		}
	}
	uint16_t nQueues = static_cast<uint16_t>(round(pow(probeRadius * resolution - sqrt(2.0), 2.0)));
	if (nQueues < 1)
		nQueues = 1;
	HierarchicalQueue* HQ1 = new HierarchicalQueue(nQueues);
	HierarchicalQueue* HQ2 = new HierarchicalQueue(nQueues);

	/* for all voxels */
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				/* create the first shell by pushing the surface voxels in the list */
				if (surface->getVoxel(i, j, k)) {
					HQ1->push(voxel(i, j, k), 0);
					/* a surface voxel has zero distance from itself */
					nbv[i][j][k] = voxel(i, j, k);
					dMap[i][j][k] = 0;
				}
				else {
					dMap[i][j][k] = UINT16_MAX;
				}
			}
		}
	}

	while (!HQ1->empty()) {
		/*current voxel*/
		voxel cv = HQ1->front();
		HQ1->pop();
		voxel curr_nbv = nbv[cv.ix][cv.iy][cv.iz];
		uint16_t sd = dMap[cv.ix][cv.iy][cv.iz];
		voxel nb; // neighbour
		bool isEnd = true;
		for (int i = 0; i < 26; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (nb.ix >= 0 && nb.ix < plength && cpkModel->getVoxel(nb)) {
				uint16_t new_sd = nb.sDistance_to(curr_nbv);
				if (new_sd < dMap[nb.ix][nb.iy][nb.iz] && new_sd < nQueues) {
					dMap[nb.ix][nb.iy][nb.iz] = new_sd;
					nbv[nb.ix][nb.iy][nb.iz] = curr_nbv;
					HQ1->push(nb, new_sd);
					isEnd = false;
				}
			}
		}
		if (isEnd && sd >= 24) {
			HQ2->push(cv, sd);
		}
	}
	while (!HQ2->empty()) {
		/*current voxel*/
		voxel cv = HQ2->front();
		HQ2->pop();
		voxel curr_nbv = nbv[cv.ix][cv.iy][cv.iz];
		// neighbor
		voxel nb;
		for (int i = 0; i < 124; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (nb.ix >= 0 && nb.ix < plength && cpkModel->getVoxel(nb)) {
				uint16_t new_sd = nb.sDistance_to(curr_nbv);
				if (new_sd < dMap[nb.ix][nb.iy][nb.iz] && new_sd < nQueues) {
					dMap[nb.ix][nb.iy][nb.iz] = new_sd;
					nbv[nb.ix][nb.iy][nb.iz] = curr_nbv;
					HQ2->push(nb, new_sd);
				}
			}
		}
	}

	/* update the volumetric model */
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (dMap[i][j][k] <= nQueues) {
					cpkModel->clearVoxel(i, j, k);
				}
			}
		}
	}

	/* Delete the partial Euclidean Distance Map. */
	for (int i = 0; i < plength; i++) {
		for (int j = 0; j < pwidth; j++) {
			delete[] nbv[i][j];
			delete[] dMap[i][j];
		}
	}
	for (int i = 0; i < plength; i++) {
		delete[] nbv[i];
		delete[] dMap[i];
	}
	delete[] nbv;
	nbv = NULL;
	delete[] dMap;
	dMap = NULL;

	delete HQ1;
	HQ1 = NULL;
	delete HQ2;
	HQ2 = NULL;

	std::cout << "Re-building surface\n";
	buildSurface_old();
} /* fastRegionGrowingEDT() */

/** Prints the 3D voxelized representation to file using the PCD
 * (Point Cloud Data) file format.
 *
 * Each PCD file contains a header that identifies and declares
 * certain properties of the point cloud data stored in the file.
 * The header of a PCD must be encoded in ASCII.
 * Storing point cloud data in both a simple ascii form with each
 * point on a line, space or tab separated, without any other
 * characters on it, as well as in a binary dump format, allows
 * us to have the best of both worlds: simplicity and speed,
 * depending on the underlying application. The ascii format
 * allows users to open up point cloud files and plot them using
 * standard software tools like gnuplot or manipulate them using
 * tools like sed, awk, etc.
 *
 * For a detailed description of the PCD (Point Cloud Data) file
 * format specification see:
 * http://pointclouds.org/documentation/tutorials/pcd_file_format.php
 *
 * \param filename	Name of the output file. The '.pcd' extension
 * 					is added automatically.
 *
 * \throws ofstream::failure
 */
void MolecularSurface::outputSurfacePCDModel(std::string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + filename + ".pcd");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream  << std::setprecision(11);
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (surface->getVoxel(i, j, k))
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}
	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n" << "FIELDS x y z\n" << "SIZE 4 4 4\n"
			<< "TYPE F F F\n" << "COUNT 1 1 1\n" << "WIDTH "
			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
			<< "\n" << "DATA ascii";

	while (!surfaceVoxels.empty()) {
		file_stream << "\n" << surfaceVoxels.front().ix / resolution - ptran.x
				<< " " << surfaceVoxels.front().iy / resolution - ptran.y
				<< " " << surfaceVoxels.front().iz / resolution - ptran.z;
		surfaceVoxels.pop();
	}
	file_stream.close();
} /* outputSurfacePCDModel() */

/**
 * Method that calculates the hydrophobicity map for the given molecule.
 * Each atom is assigned the hydrophobicity value of the residue it belongs.
 * The algorithm is very similar to the volumetric model creation process:
 * the only difference is that instead of setting voxels to 1 if they are
 * covered by some atom, we add the atom's hydrophobicity value to all voxels
 * of me molecular surface it covers.
 *
 * @param hydrophobicity	the hydrophobicity scale to use in the calculation
 */
void MolecularSurface::calculateHydrophobicity(std::map<std::string, float> const & hydrophobicity) {
	/* Initialize the hydrophobicity map*/
	hydrophobicityMap.resize(plength);
	for (int i = 0; i < plength; ++i) {
		hydrophobicityMap[i].resize(pwidth);
		for (int j = 0; j < pwidth; ++j) {
			hydrophobicityMap[i][j].resize(pheight);
			for (int k = 0; k < pheight; ++k) {
				hydrophobicityMap[i][j][k] = 0;
			}
		}
	}
	uint16_t cx, cy, cz; /**< Discretized coordinates of the atom's center. */
	uint16_t radius; /**< discretized atomic radius */
	double hphob; /**< hydrophobicity value for the current atom */

	uint16_t y_max; /**< maximum value for the squared y coordinate */
	uint16_t z_max; /**< maximum value for the squared z coordinate */

	std::vector<atom>::const_iterator atm;
	for (atm = atoms->begin(); atm != atoms->end(); ++atm) {
		if (atm->hasZeroRadius())
			continue;
		map<string, float>::const_iterator entry = hydrophobicity.find(atm->residue_name);
		if (entry != hydrophobicity.end())
			hphob = entry->second;
		else
			continue;
		/* Translate and discretize the coordinates */

		cx = static_cast<uint16_t>((atm->x + ptran.x) * resolution + 0.5);
		cy = static_cast<uint16_t>((atm->y + ptran.y) * resolution + 0.5);
		cz = static_cast<uint16_t>((atm->z + ptran.z) * resolution + 0.5);

		radius = static_cast<uint16_t>((atm->radius + probeRadius) * resolution + 0.5);

		// x^2 + y^2 + z^2 = r^2
		for (int x = -radius; x <= radius; ++x) {
			// y^2 + z^2 = r^2 - x^2
			y_max = static_cast<int>(sqrt(radius * radius - x * x) + 0.5);
			for (int y = -y_max; y <= y_max; ++y) {
				// z^2 = r^2 - x^2 - y^2
				z_max = static_cast<int>(sqrt(radius * radius - x * x - y * y) + 0.5);
				for (int z = -z_max; z <= z_max; ++z) {
					if (surface->getVoxel(cx + x, cy + y, cz + z))
						hydrophobicityMap[cx + x][cy + y][cz + z] += hphob
								/ sqrt(x * x + y * y + z * z);
				}
			}
		}
	}
	/* Now we normalize the hydrophobicity values between 0 and 1 */
	normalizePositive(hydrophobicityMap);
	normalizeNegative(hydrophobicityMap);

} /* calculateHydrophobicity() */
void MolecularSurface::calculateHydrophobicity(std::map<std::string, float> const & hydrophobicity, size_t window_size) {
	/* Initialize the hydrophobicity map*/
	hydrophobicityMap.resize(plength);
	for (int i = 0; i < plength; ++i) {
		hydrophobicityMap[i].resize(pwidth);
		for (int j = 0; j < pwidth; ++j) {
			hydrophobicityMap[i][j].resize(pheight);
			for (int k = 0; k < pheight; ++k) {
				hydrophobicityMap[i][j][k] = 0;
			}
		}
	}
	vector<pair<atom,double>> atoms_hphob;
	for (auto const & atm : *atoms) {
		atoms_hphob.push_back(make_pair(atm, 0.0));
	}
	sort(atoms_hphob.begin(), atoms_hphob.end(),
			[](pair<atom,double> const & l, pair<atom,double> const & r) {return l.first.residue_number < r.first.residue_number;});
	for (size_t ii = 0; ii < atoms_hphob.size(); ++ii) {
		double average_hphob;
		for (size_t jj = ii - (window_size - 1) / 2; jj <= ii + (window_size - 1)/2; ++jj) {
			if (jj < 0 || jj >= atoms_hphob.size())
				continue;
			auto entry = hydrophobicity.find(atoms_hphob[jj].first.residue_name);
				if (entry != hydrophobicity.end())
					average_hphob += entry->second;
		}
		average_hphob /= window_size;
		atoms_hphob[ii].second = average_hphob;
	}
	uint16_t cx, cy, cz; /**< Discretized coordinates of the atom's center. */
	uint16_t radius; /**< discretized atomic radius */

	uint16_t y_max; /**< maximum value for the squared y coordinate */
	uint16_t z_max; /**< maximum value for the squared z coordinate */

	for (auto const & atm_h : atoms_hphob) {
		if (atm_h.first.hasZeroRadius())
			continue;
		/* Translate and discretize the coordinates */
		cx = static_cast<uint16_t>((atm_h.first.x + ptran.x) * resolution + 0.5);
		cy = static_cast<uint16_t>((atm_h.first.y + ptran.y) * resolution + 0.5);
		cz = static_cast<uint16_t>((atm_h.first.z + ptran.z) * resolution + 0.5);

		radius = static_cast<uint16_t>((atm_h.first.radius + probeRadius) * resolution + 0.5);

		// x^2 + y^2 + z^2 = r^2
		for (int x = -radius; x <= radius; ++x) {
			// y^2 + z^2 = r^2 - x^2
			y_max = static_cast<int>(sqrt(radius * radius - x * x) + 0.5);
			for (int y = -y_max; y <= y_max; ++y) {
				// z^2 = r^2 - x^2 - y^2
				z_max = static_cast<int>(sqrt(radius * radius - x * x - y * y) + 0.5);
				for (int z = -z_max; z <= z_max; ++z) {
					if (surface->getVoxel(cx + x, cy + y, cz + z))
						hydrophobicityMap[cx + x][cy + y][cz + z] += atm_h.second
								/ sqrt(x * x + y * y + z * z);
				}
			}
		}
	}
	/* Now we normalize the hydrophobicity values between 0 and 1 */
//	normalizePositive(hydrophobicityMap);
//	normalizeNegative(hydrophobicityMap);

} /* calculateHydrophobicity() */
void MolecularSurface::calculateAtomHydrophobicity(std::unordered_map<pair<string, string>, int> const & hydrophobicity) {
	/* Initialize the hydrophobicity map*/
	hydrophobicityMap.resize(plength);
	for (int i = 0; i < plength; ++i) {
		hydrophobicityMap[i].resize(pwidth);
		for (int j = 0; j < pwidth; ++j) {
			hydrophobicityMap[i][j].resize(pheight);
			for (int k = 0; k < pheight; ++k) {
				hydrophobicityMap[i][j][k] = 0;
			}
		}
	}
	uint16_t cx, cy, cz; /**< Discretized coordinates of the atom's center. */
	uint16_t radius; /**< discretized atomic radius */
	int hphob; /**< hydrophobicity value for the current atom */

	uint16_t y_max; /**< maximum value for the squared y coordinate */
	uint16_t z_max; /**< maximum value for the squared z coordinate */

	std::vector<atom>::const_iterator atm;
	for (atm = atoms->begin(); atm != atoms->end(); ++atm) {
		if (atm->hasZeroRadius())
			continue;
		auto entry = hydrophobicity.find(make_pair(atm->residue_name, trim(atm->atom_name)));
		if (entry != hydrophobicity.end())
			hphob = -1 *  entry->second;
		else
			continue;
		/* Translate and discretize the coordinates */
		cx = static_cast<uint16_t>((atm->x + ptran.x) * resolution + 0.5);
		cy = static_cast<uint16_t>((atm->y + ptran.y) * resolution + 0.5);
		cz = static_cast<uint16_t>((atm->z + ptran.z) * resolution + 0.5);

		radius = static_cast<uint16_t>((atm->radius + probeRadius) * resolution + 0.5);

		// x^2 + y^2 + z^2 = r^2
		for (int x = -radius; x <= radius; ++x) {
			// y^2 + z^2 = r^2 - x^2
			y_max = static_cast<int>(sqrt(radius * radius - x * x) + 0.5);
			for (int y = -y_max; y <= y_max; ++y) {
				// z^2 = r^2 - x^2 - y^2
				z_max = static_cast<int>(sqrt(radius * radius - x * x - y * y) + 0.5);
				for (int z = -z_max; z <= z_max; ++z) {
					if (surface->getVoxel(cx + x, cy + y, cz + z) && !(x == 0 && y == 0 && z == 0))
						hydrophobicityMap[cx + x][cy + y][cz + z] += hphob
								/ (float) sqrt(x * x + y * y + z * z);
				}
			}
		}
	}
	/* Now we normalize the hydrophobicity values between 0 and 1 */
//	normalizePositive(hydrophobicityMap);
//	normalizeNegative(hydrophobicityMap);
}




/**
 * This method outputs the molecular surface with the hydrophobicity
 * attribute mapped onto it in the Point Cloud Data format.
 *
 * The hydrophobicity is represented in terms of colors of the surface.
 * The blue color indicates a low hydrophobicity (hydrophily) of the
 * region, while the red color indicates a high level of hydrophobicity.
 *
 * \param filename	Name of the output file. The '-hydrophobicity.pcd'
 * 					suffix is added automatically.
 * \throws ofstream::failure
 */
void MolecularSurface::outputHydrophobicity(std::string const & filename) {
	assert(!hydrophobicityMap.empty());
	assert(surface != NULL);

	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + filename + "-hydrophobicity.pcd");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream  << std::setprecision(11);
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (surface->getVoxel(i, j, k))
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}

	union {
		uint8_t rgb_byte[4];
		float rgb_flt;
	} color_RGB;

	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n" << "FIELDS x y z rgb\n" << "SIZE 4 4 4 "<< sizeof(color_RGB) <<"\n"
			<< "TYPE F F F F\n" << "COUNT 1 1 1 1\n" << "WIDTH "
			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
			<< "\n" << "DATA ascii";
	float hphob;

	int x, y, z;
	uint8_t r, g, b;

	while (!surfaceVoxels.empty()) {
		x = surfaceVoxels.front().ix;
		y = surfaceVoxels.front().iy;
		z = surfaceVoxels.front().iz;

//		hphob = (hydrophobicityMap[x][y][z] - 0.5) * 2;
		hphob = hydrophobicityMap[x][y][z];

		if (hphob <= 0) {
			r = 255;
			g = (uint8_t) floor((1.0 + hphob) * 255);
			b = (uint8_t) floor((1.0 + hphob) * 255);
		} else {
			r = (uint8_t) floor((1.0 - hphob) * 255);
			g = (uint8_t) floor((1.0 - hphob) * 255);
			b = 255;
		}

		color_RGB.rgb_byte[0] = b;
		color_RGB.rgb_byte[1] = g;
		color_RGB.rgb_byte[2] = r;
		color_RGB.rgb_byte[3] = 0;

		file_stream << "\n" << x / resolution - ptran.x << " "
				<< y / resolution - ptran.y << " "
				<< z / resolution - ptran.z << " "
				<< color_RGB.rgb_flt;
		surfaceVoxels.pop();
	}
	file_stream.close();
} /* outputHydrophobicity() */

/**
 * This method calculates electrostatic potential according to
 * Coulomb's law:
 *
 * 						φ = Σ [qi / (εdi)]
 *
 * φ is the potential (which varies in space), q are the atomic
 * partial charges, d are the distances from the atoms, and ε is
 * the dielectric, representing screening by the medium or solvent.
 * A distance-dependent dielectric (ε = Cd where C is some constant)
 * is sometimes used to approximate screening by implicit solvent.
 */
void MolecularSurface::calculateCoulombicPotentials() {
	/* Initialize the coulombicPotentials map*/
	coulombicPotentials.resize(plength);
	for (int i = 0; i < plength; ++i) {
		coulombicPotentials[i].resize(pwidth);
		for (int j = 0; j < pwidth; ++j) {
			coulombicPotentials[i][j].resize(pheight);
			for (int k = 0; k < pheight; ++k) {
				coulombicPotentials[i][j][k] = 0;
			}
		}
	}

	float currentPotential, currentDistance;
	float cutoff = 6; // 6 Ångstrom cutoff value
	uint16_t radius = static_cast<uint16_t>(cutoff * resolution + 0.5);
	uint16_t y_max; /** maximum value for the squared y coordinate */
	uint16_t z_max; /**< maximum value for the squared z coordinate */
	int cx, cy, cz; /**< Current atom center coordinates. */
	int si, sj, sk; /**< Surrounding voxels' coordinates. */
	std::vector<atom>::const_iterator atm;
	for (atm = atoms->begin(); atm != atoms->end(); ++atm) {
		currentPotential = atm->charge;
		cx = static_cast<uint16_t>((atm->x + ptran.x) * resolution + 0.5);
		cy = static_cast<uint16_t>((atm->y + ptran.y) * resolution + 0.5);
		cz = static_cast<uint16_t>((atm->z + ptran.z) * resolution + 0.5);

		// x^2 + y^2 + z^2 = r^2
		for (int x = -radius; x <= radius; ++x) {
			// y^2 + z^2 = r^2 - x^2
			y_max = static_cast<int>(sqrt(radius * radius - x * x) + 0.5);
			for (int y = -y_max; y <= y_max; ++y) {
				// z^2 = r^2 - x^2 - y^2
				z_max = static_cast<int>(sqrt(radius * radius - x * x - y * y) + 0.5);
				for (int z = -z_max; z <= z_max; ++z) {
					si = cx + x;
					sj = cy + y;
					sk = cz + z;
					if (sk > -1 && sk < pheight && sj > -1 && sj < pwidth && si > -1 && si < plength) {
						if (x == 0 && y == 0 && z == 0)
							coulombicPotentials[si][sj][sk] += currentPotential;
						else {
							currentDistance = sqrt(x * x + y * y + z * z) / resolution;
							coulombicPotentials[si][sj][sk] += currentPotential / (EPSILON * currentDistance);
						}

					}
				}
			}
		}
	}
	/* Now we normalize the coulombic potentials values between 0 and 1 */
	normalizeGrid(coulombicPotentials);
} /* calculateCoulombicPotentials() */

/**
 * This method outputs the molecular surface with the Coulombic potential
 * attribute mapped onto it in the Point Cloud Data format.
 *
 * The Coulombic potential is represented in terms of colors of the surface.
 * The blue color indicates a low electrostatic potential of the
 * region, while the red color indicates a high Coulombic potential.
 *
 * \param filename	Name of the output file. The '-coulombicPotentials.pcd'
 * 					suffix is added automatically.
 * \throws ofstream::failure
 */
void MolecularSurface::outputCoulombicPotentials(std::string const & filename) {
	assert(!coulombicPotentials.empty());
	assert(surface != NULL);
	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + filename + "-coulombicPotentials.pcd");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream  << std::setprecision(11);
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (surface->getVoxel(i, j, k)) {
					surfaceVoxels.push(voxel(i, j, k));
				}
			}
		}
	}

	union {
		uint8_t rgb_byte[4];
		float rgb_flt;
	} color_RGB;

	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n" << "FIELDS x y z rgb\n" << "SIZE 4 4 4 "<< sizeof(color_RGB) <<"\n"
			<< "TYPE F F F F\n" << "COUNT 1 1 1 1\n" << "WIDTH "
			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
			<< "\n" << "DATA ascii";
	float cp; /**< current coulombic potential */
	int x, y, z;
	int r, g, b;
	while (!surfaceVoxels.empty()) {
		x = surfaceVoxels.front().ix;
		y = surfaceVoxels.front().iy;
		z = surfaceVoxels.front().iz;
		cp = coulombicPotentials[x][y][z];

		if (cp <= 0) {
			r = 255;
			g = (uint8_t) floor((1.0 + cp) * 255);
			b = (uint8_t) floor((1.0 + cp) * 255);
		} else {
			r = (uint8_t) floor((1.0 - cp) * 255);
			g = (uint8_t) floor((1.0 - cp) * 255);
			b = 255;
		}

		color_RGB.rgb_byte[0] = b;
		color_RGB.rgb_byte[1] = g;
		color_RGB.rgb_byte[2] = r;
		color_RGB.rgb_byte[3] = 0;

		file_stream << "\n" << x / resolution - ptran.x << " "
				<< y / resolution - ptran.y << " "
				<< z / resolution - ptran.z << " "
				<< color_RGB.rgb_flt;
		surfaceVoxels.pop();
	}
	file_stream.close();
} /* outputCoulombicPotentials() */

/**
 * This method creates the DX potentials grid from reading an APBS-generated
 * input openDX file. The positive and negative potentials are stored in two
 * different grids, enabling us to use the Zernike invariants during the
 * docking process.
 *
 * @param input_DX		input openDX filename
 */
void MolecularSurface::calculateDXPotentials(const std::string & input_DX) {
	assert(surface != NULL);
	DXPotentials = new PotentialGridDX(input_DX.c_str());
	size_t gridSize = plength * pwidth * pheight;
	potentialsDX = array3D(plength, array2D(pwidth, array1D(pheight, 0.0)));
//	negativeDXPotentials = array3D(plength, array2D(pwidth, array1D(pheight, 0.0)));

	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (surface->getVoxel(i, j, k)) {
					point3D p(i, j, k);
					p /= resolution;
					p -= ptran;
					double pDX = DXPotentials->getPotential(p);
					potentialsDX[i][j][k] = pDX;

//					if (pDX < 0)
//						negativeDXPotentials[i][j][k] = -pDX;
//					else
//						positiveDXPotentials[i][j][k] = pDX;
				}
			}
		}
	}
//	normalizePositive(positiveDXPotentials);
//	normalizePositive(negativeDXPotentials);
} /* calculateDXPotentials() */
/**
 * This method outputs the molecular surface with the APBS-calculated
 * potential attribute mapped onto it in the Point Cloud Data format.
 *
 * \param filename	Name of the output file. The '-DXPotentials.pcd'
 * 					suffix is added automatically.
 * \throws ofstream::failure
 */
void MolecularSurface::outputDXPotentials(std::string const & outname) {
	assert(surface != NULL);
	assert(!potentialsDX.empty());

	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + outname + "-DXPotentials.pcd");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream  << std::setprecision(11);
	std::queue<voxel> surfacePoints;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (surface->getVoxel(i, j, k)) {
					surfacePoints.push(voxel(i, j, k));
				}
			}
		}
	}

	union {
		uint8_t rgb_byte[4];
		float rgb_flt;
	} color_RGB;

	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n" << "FIELDS x y z rgb\n" << "SIZE 4 4 4 "<< sizeof(color_RGB) <<"\n"
			<< "TYPE F F F F\n" << "COUNT 1 1 1 1\n" << "WIDTH "
			<< surfacePoints.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfacePoints.size()
			<< "\n" << "DATA ascii";
	float pDX; /**< current potential */
	int x, y, z;
	uint8_t r, g, b;

	while (!surfacePoints.empty()) {
		x = surfacePoints.front().ix;
		y = surfacePoints.front().iy;
		z = surfacePoints.front().iz;

		pDX = potentialsDX[x][y][z];//positiveDXPotentials[x][y][z] - negativeDXPotentials[x][y][z];
		if (pDX <= 0) {
			r = 255;
			g = (uint8_t) floor((1.0 + pDX) * 255);
			b = (uint8_t) floor((1.0 + pDX) * 255);
		} else {
			r = (uint8_t) floor((1.0 - pDX) * 255);
			g = (uint8_t) floor((1.0 - pDX) * 255);
			b = 255;
		}

		color_RGB.rgb_byte[0] = b;
		color_RGB.rgb_byte[1] = g;
		color_RGB.rgb_byte[2] = r;
		color_RGB.rgb_byte[3] = 0;

		file_stream << "\n" << x / resolution - ptran.x << " "
				<< y / resolution - ptran.y << " "
				<< z / resolution - ptran.z << " "
				<< color_RGB.rgb_flt;
		surfacePoints.pop();
	}
	file_stream.close();
} /* outputDXPotentials() */
/**
 * Very simple method for extracting the patch centers from the current molecular
 * surface. In fact it is a very simple decimation algorithm. Scanning the 3D
 * grid, the first surface point encountered becomes the first patch center.
 * At this voxel, a sphere of 'minDistance' radius is constructed, and all the
 * voxels inside this sphere are 'cleared' from the search. The scanning restarts,
 * until a new voxel is encountered (which has not been cleared yet). This point
 * will be marked as the second patch center, a new sphere with 'minDistance' radius
 * centered on this voxel will be created and all the voxels inside the sphere will
 * be cleared. The algorithm proceeds until the whole 3D grid containing the molecular
 * surface has been scanned and there are no 'uncleared' voxels left.
 *
 * @param minDistance	minimum distance between two different patch centers (in Å).
 */
void MolecularSurface::extractPatchCenters(float minDistance) {
	assert(surface != NULL);
	/* discretized minimum distance between two different surface patches */
	float minD = minDistance * resolution;
	float s_minD = minD * minD;
	uint16_t d_minDistance = (uint16_t) ceil(minD);


	std::list<voxel_offset> sphereOffsets;
	get_sphere_offsets(d_minDistance, sphereOffsets);

	voxelGrid * temp = new voxelGrid(*surface);
	if (!patchCenters.empty())
		patchCenters.clear();
	int x, y, z;

	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (temp->getVoxel(i, j, k)) { // && is_SA(voxel(i, j, k))){// && hydrophobicityMap[i][j][k] > 0) {
					/* put the current voxel in the patchCenters list */
					patchCenters.push_back(voxel(i, j, k));
					/* clear all surrounding voxels within minDistance from voxel (i, j, k) */
					for (auto const & o : sphereOffsets) {
						x = i + o.i;
						y = j + o.j;
						z = k + o.k;
						if (x > -1 && x < plength && y > -1 && y < pwidth && z > -1 && z < pheight) {
							temp->clearVoxel(x, y, z);
						}
					}
				} // if
			} // for k
		} // for j
	} // for i
	delete temp;
} /* extractSurfacePatches() */


void MolecularSurface::extractPatchCenters(float minDistance, vector<voxel> & patchCenters) {
	assert(surface != NULL);
	/* discretized minimum distance between two different surface patches */
	float minD = minDistance * resolution;
	float s_minD = minD * minD;
	uint16_t d_minDistance = (uint16_t) ceil(minD);


	std::list<voxel_offset> sphereOffsets;
	get_sphere_offsets(d_minDistance, sphereOffsets);

	voxelGrid * temp = new voxelGrid(*surface);
	int x, y, z;

	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (temp->getVoxel(i, j, k)) { // && is_SA(voxel(i, j, k))){// && hydrophobicityMap[i][j][k] > 0) {
					/* put the current voxel in the patchCenters list */
					patchCenters.push_back(voxel(i, j, k));
					/* clear all surrounding voxels within minDistance from voxel (i, j, k) */
					for (auto const & o : sphereOffsets) {
						x = i + o.i;
						y = j + o.j;
						z = k + o.k;
						if (x > -1 && x < plength && y > -1 && y < pwidth && z > -1 && z < pheight) {
							temp->clearVoxel(x, y, z);
						}
					}
				} // if
			} // for k
		} // for j
	} // for i
	delete temp;
} /* extractSurfacePatches() */

void MolecularSurface::extractInterfacePatchCenters(float minDistance, array3D const & interface) {
	assert(surface != NULL);
	/* discretized minimum distance between two different surface patches */
	float minD = minDistance * resolution;
	float s_minD = minD * minD;
	uint16_t d_minDistance = (uint16_t) ceil(minD);


	std::list<voxel_offset> sphereOffsets;
	get_sphere_offsets(d_minDistance, sphereOffsets);

	voxelGrid * temp = new voxelGrid(*surface);
	if (!patchCenters.empty())
		patchCenters.clear();
	int x, y, z;

	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (temp->getVoxel(i, j, k)) {
					/* put the current voxel in the patchCenters list */
					if (interface[i][j][k] > 0)
						patchCenters.push_back(voxel(i, j, k));
					/* clear all surrounding voxels within minDistance from voxel (i, j, k) */
					for (auto const & o : sphereOffsets) {
						x = i + o.i;
						y = j + o.j;
						z = k + o.k;
						if (x > -1 && x < plength && y > -1 && y < pwidth && z > -1 && z < pheight) {
							temp->clearVoxel(x, y, z);
						}
					}
				} // if
			} // for k
		} // for j
	} // for i
	delete temp;
} /* extractSurfacePatches() */
void MolecularSurface::extractInterfacePatchCenters(float minDistance, array3D const & interface, vector<voxel> & patchCenters) {
	assert(surface != NULL);
	/* discretized minimum distance between two different surface patches */
	float minD = minDistance * resolution;
	float s_minD = minD * minD;
	uint16_t d_minDistance = (uint16_t) ceil(minD);


	std::list<voxel_offset> sphereOffsets;
	get_sphere_offsets(d_minDistance, sphereOffsets);

	voxelGrid * temp = new voxelGrid(*surface);
	if (!patchCenters.empty())
		patchCenters.clear();
	int x, y, z;

	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (temp->getVoxel(i, j, k)) {
					/* put the current voxel in the patchCenters list */
					if (interface[i][j][k] > 0)
						patchCenters.push_back(voxel(i, j, k));
					/* clear all surrounding voxels within minDistance from voxel (i, j, k) */
					for (auto const & o : sphereOffsets) {
						x = i + o.i;
						y = j + o.j;
						z = k + o.k;
						if (x > -1 && x < plength && y > -1 && y < pwidth && z > -1 && z < pheight) {
							temp->clearVoxel(x, y, z);
						}
					}
				} // if
			} // for k
		} // for j
	} // for i
	delete temp;
} /* extractSurfacePatches() */
/**
 * Debugging method.
 * Prints the extracted patch centers to PCD.
 */
void MolecularSurface::outputPatchCentersPCDModel(std::string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open("./output/" + filename + "-patch_centers.pcd");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream  << std::setprecision(11);
	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n"
			<< "FIELDS x y z\n"
			<< "SIZE 4 4 4\n" << "TYPE F F F\n"
			<< "COUNT 1 1 1\n" << "WIDTH " << patchCenters.size()
			<< "\n" << "HEIGHT 1\n" << "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS "
			<< patchCenters.size() << "\n" << "DATA ascii";
	for (auto const & center : patchCenters) {
		file_stream << "\n"
				<< center.ix / resolution - ptran.x << " "
				<< center.iy / resolution - ptran.y << " "
				<< center.iz / resolution - ptran.z;
	}
	file_stream.close();
} /* outputPatchCentersPCDModel() */
/**
 * Debugging method.
 * Prints the extracted patch centers to PCD.
 */
//void MolecularSurface::outputPatchCentersNormalsPCDModel(std::string const & filename, float patchRadius, float normalRadius) {
//	std::ofstream file_stream;
//	makeDirectory("./output");
//	file_stream.open("./output/" + filename + "-patch_centers.pcd");
//	if (!file_stream.is_open()) {
//		throw std::ofstream::failure("Error opening output file.");
//	}
//	file_stream  << std::setprecision(11);
//	/* File header */
//	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
//			<< "VERSION .7\n"
//			<< "FIELDS x y z normal_x normal_y normal_z curvature\n"
//			<< "SIZE 4 4 4 4 4 4 4\n" << "TYPE F F F F F F F\n"
//			<< "COUNT 1 1 1 1 1 1 1\n" << "WIDTH " << patchCenters.size()
//			<< "\n" << "HEIGHT 1\n" << "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS "
//			<< patchCenters.size() << "\n" << "DATA ascii";
//	uint16_t counter = 0;
//	for (auto const & center : patchCenters) {
//		SurfacePatch currentPatch(counter, patchRadius, normalRadius, center, resolution, surface, cpkModel, ptran);
//		file_stream << "\n"
//				<< currentPatch.patchCenter.ix / resolution - ptran.x << " "
//				<< currentPatch.patchCenter.iy / resolution - ptran.y << " "
//				<< currentPatch.patchCenter.iz / resolution - ptran.z << " "
//				<< currentPatch.patchNormal.x << " "
//				<< currentPatch.patchNormal.y << " "
//				<< currentPatch.patchNormal.z << " ";
//		++counter;
//	}
//	file_stream.close();
//} /* outputPatchCentersPCDModel() */
void MolecularSurface::calculateSurfaceDescriptors(float patchRadius, int maxOrder,
		array3D const & interface, std::vector<CompactPatchDescriptor> & descriptors) {
	if (!descriptors.empty())
		descriptors.clear();
	size_t n = patchCenters.size();
	descriptors.resize(n);
	/* the id corresponds to the index of the current descriptor in the output vector */
#pragma omp parallel for
	for (size_t i = 0; i < n; ++i) {
		SurfacePatch currentPatch(i, patchRadius, patchCenters[i], resolution,
				surface, cpkModel, potentialsDX, hydrophobicityMap, ptran,
				interface, BLAM930101_, BIOV880101_, MAXF760101_,
				TSAJ990101_, NAKH920108_, CEDJ970104_, LIFS790101_,
				MIYS990104_);
		currentPatch.calculateZernikeDescriptors(maxOrder);
		descriptors[i] = CompactPatchDescriptor(currentPatch);
	}
}


void MolecularSurface::calculateSurfaceDescriptors(float patchRadius, int maxOrder, vector<voxel> const & patchCenters,
		array3D const & interface, std::vector<CompactPatchDescriptor> & descriptors) {
	if (!descriptors.empty())
		descriptors.clear();
	size_t n = patchCenters.size();
	descriptors.resize(n);
	/* the id corresponds to the index of the current descriptor in the output vector */
#pragma omp parallel for
	for (size_t i = 0; i < n; ++i) {
		SurfacePatch currentPatch(i, patchRadius, patchCenters[i], resolution,
				surface, cpkModel, potentialsDX, hydrophobicityMap, ptran,
				interface, BLAM930101_, BIOV880101_, MAXF760101_,
				TSAJ990101_, NAKH920108_, CEDJ970104_, LIFS790101_,
				MIYS990104_);
		currentPatch.calculateZernikeDescriptors(maxOrder);
		descriptors[i] = CompactPatchDescriptor(currentPatch);
	}
}
/**
 * Calculate the accessible-surface area of atoms.
 * Uses the simple Shrake-Rupley algorithm, that generates a
 * relatively uniform density of dots over every atoms and
 * eliminates those within the sphere of another atom. The remaining
 * dots are used to calculate the area.
 *
 * Reference: A. Shrake & J. A. Rupley. "Environment and Exposure to
 * Solvent of Protein Atoms. Lysozyme and Insulin." J Mol Biol. 79
 * (1973) 351- 371.
 * @param n_sphere_points	the number of sphere points to use for the
 * 							ASA calculation
 * @param per_atom_asa		vector containing the ASA of each atom; element i
 * 							gives the area of
 * @param per_atom_sap		vector containing the surface accessible points for
 * 							each atom
 *
 * @return	the Accessible Surface Area of the current molecule
 */
double MolecularSurface::calculateAccessibleSurfaceArea(std::vector<double> & per_atom_asa,
		std::vector<std::vector<point3D>> & per_atom_sap, size_t n_sphere_points) {
	assert(cpkModel != NULL);
	size_t num_atoms = atoms->size();
	if (num_atoms == 0)
		return 0;
	per_atom_asa.resize(num_atoms);
	per_atom_sap.resize(num_atoms);


	if (sphere_points.empty())
		generateSpherePoints(n_sphere_points);

	double total_ASA = 0;

	double c = 4.0 * M_PI / n_sphere_points;
	for (size_t i = 0; i < num_atoms; ++i) {
		if ((*atoms)[i].hasZeroRadius())
			per_atom_asa[i] = 0;
		else {
			double radius = (*atoms)[i].radius + probeRadius;
			size_t n_accessible_pts = 0;
			point3D atm_center((*atoms)[i]);
			for (size_t j = 0; j < n_sphere_points; ++j) {
				point3D p = radius * sphere_points[j] + atm_center;
				point3D scaled_p = resolution * (p + ptran);
				uint16_t d_px = static_cast<uint16_t>(scaled_p.x + 0.5);
				uint16_t d_py = static_cast<uint16_t>(scaled_p.y + 0.5);
				uint16_t d_pz = static_cast<uint16_t>(scaled_p.z + 0.5);
				if (is_SA(voxel(d_px, d_py, d_pz))) {
					++n_accessible_pts;
					per_atom_sap[i].push_back(p);
				}
			}
			per_atom_asa[i] = c * n_accessible_pts * radius * radius;
			total_ASA += per_atom_asa[i];
		}
	}
	return total_ASA;
}
/**
 * Solvation energy of a complex is strongly correlated to its binding free energy.
 * ASP (Atom Solvation Parameters) model, first proposed by Eisenberg and McLachlan
 * in 1986 (Eisenberg D, McLachlan A: Solvation energy in protein folding and binding.
 * Nature 1986.), is one of the most successful models for solvation energy calculation.
 * ASP model assumes that the solvation energy of an atom or an atom-group is
 * proportional to the area of its solvent accessible surface and so the total solvation
 * energy of a molecule is:
 *
 * 		ΔG = ∑ σi*Ai
 *
 * where Ai is the solvent-accessible surface area (ASA) of atom i and σi is the ASP
 * value of atom i, which can be determined experimentally.
 *
 * @param asp			reference to the map containing the atom solvation parameters
 * 						to use for the solvation energy calculation
 * @param per_atom_asa 	vector containing the ASA of each atom in the molecule
 * @return	the solvation energy of the current molecule
 */
double MolecularSurface::calculateSolvationEnergy(std::map<std::string, float> const & asp,
		std::vector<double> const & per_atom_asa) {
	double dG = 0;
	size_t num_atoms = atoms->size();
	if (num_atoms == 0)
		return 0;
	assert (per_atom_asa.size() == num_atoms);

	float charge;
	std::string atom_name;
	float sigma;

	for (size_t i = 0; i < num_atoms; ++i) {
		charge = (*atoms)[i].charge;
		atom_name = trim((*atoms)[i].atom_name);

		if (atom_name[0] == 'C')
			sigma = asp.find("C")->second;
		else if (atom_name[0] == 'O')
			if (charge < 0)
				sigma = asp.find("O-")->second;
			else
				sigma = asp.find("O/N")->second;
		else if (atom_name[0] == 'N')
			if (charge > 0)
				sigma = asp.find("N+")->second;
			else
				sigma = asp.find("O/N")->second;
		else if (atom_name[0] == 'S')
			sigma = asp.find("S")->second;
		else
			sigma = 0;

		dG += per_atom_asa[i] * sigma;
	}
	return dG;
}
/**
 * Solvation energy of a complex is strongly correlated to its binding free energy.
 * ASP (Atom Solvation Parameters) model, first proposed by Eisenberg and McLachlan
 * in 1986 (Eisenberg D, McLachlan A: Solvation energy in protein folding and binding.
 * Nature 1986.), is one of the most successful models for solvation energy calculation.
 * ASP model assumes that the solvation energy of an atom or an atom-group is
 * proportional to the area of its solvent accessible surface and so the total solvation
 * energy of a molecule is:
 *
 * 		ΔG = ∑ σi*Ai
 *
 * where Ai is the solvent-accessible surface area (ASA) of atom i and σi is the ASP
 * value of atom i, which can be determined experimentally.
 *
 * Very similar to the previous method, the only difference is that this method calls
 * the ASA calculation procedure inside its code.
 *
 * @param asp			reference to the map containing the atom solvation parameters
 * 						to use for the solvation energy calculation
 * @return	the solvation energy of the current molecule
 */
double MolecularSurface::calculateSolvationEnergy(std::map<std::string, float> const & asp) {
	double dG = 0;
	size_t num_atoms = atoms->size();
	if (num_atoms == 0)
		return 0;

	std::vector<double> per_atom_asa;
	std::vector<std::vector<point3D>> per_atom_sap;

	calculateAccessibleSurfaceArea(per_atom_asa, per_atom_sap);

	assert (per_atom_asa.size() == num_atoms);

	float charge;
	std::string atom_name;
	float sigma;

	for (size_t i = 0; i < num_atoms; ++i) {
		charge = (*atoms)[i].charge;
		atom_name = (*atoms)[i].atom_name;

		if (atom_name[0] == 'C')
			sigma = asp.find("C")->second;
		else if (atom_name[0] == 'O')
			if (charge < 0)
				sigma = asp.find("O-")->second;
			else
				sigma = asp.find("O/N")->second;
		else if (atom_name[0] == 'N')
			if (charge > 0)
				sigma = asp.find("N+")->second;
			else
				sigma = asp.find("O/N")->second;
		else if (atom_name[0] == 'S')
			sigma = asp.find("S")->second;
		else
			sigma = 0;

		dG += per_atom_asa[i] * sigma;
	}
	return dG;
}

double MolecularSurface::calculateHydrophobicEnergy(std::map<std::string, float> const & hydrophobicity,
		std::vector<double> const & per_atom_asa) {
	double Eh = 0;
	size_t num_atoms = atoms->size();
	if (num_atoms == 0)
		return 0;
	assert (per_atom_asa.size() == num_atoms);

	std::string residue;
	float h;

	for (size_t i = 0; i < num_atoms; ++i) {
		map<string, float>::const_iterator entry = hydrophobicity.find((*atoms)[i].residue_name);
		if (entry != hydrophobicity.end())
			h = entry->second;
		else
			continue;
		Eh += per_atom_asa[i] * h;
	}
	return Eh;
}

//void MolecularSurface::outputSurfacePatchPCDModel(size_t id, float patchRadius, string const & outname) {
//	SurfacePatch currentPatch(id, patchRadius, patchCenters[id], resolution, surface, positiveDXPotentials, negativeDXPotentials, ptran);
//	if (DXPotentials != NULL)
//		currentPatch.exportPatchPotentialsPCDModel(outname);
//	else
//		currentPatch.exportPatchPCDModel(outname);
//}

void MolecularSurface::calculateHQI8() {
	BLAM930101_ = array3D(plength, array2D(pwidth, array1D(pheight, 0)));
	BIOV880101_ = array3D(plength, array2D(pwidth, array1D(pheight, 0)));
	MAXF760101_	= array3D(plength, array2D(pwidth, array1D(pheight, 0)));
	TSAJ990101_	= array3D(plength, array2D(pwidth, array1D(pheight, 0)));
	NAKH920108_	= array3D(plength, array2D(pwidth, array1D(pheight, 0)));
	CEDJ970104_	= array3D(plength, array2D(pwidth, array1D(pheight, 0)));
	LIFS790101_	= array3D(plength, array2D(pwidth, array1D(pheight, 0)));
	MIYS990104_	= array3D(plength, array2D(pwidth, array1D(pheight, 0)));

	uint16_t cx, cy, cz; /**< Discretized coordinates of the atom's center. */
	uint16_t radius; /**< discretized atomic radius */
	vector<float> hqi8; /**< hydrophobicity value for the current atom */

	uint16_t y_max; /**< maximum value for the squared y coordinate */
	uint16_t z_max; /**< maximum value for the squared z coordinate */

	std::vector<atom>::const_iterator atm;
	for (auto const & atm : *atoms) {
		if (atm.hasZeroRadius())
			continue;
		hqi8 = get_HQI8(trim(atm.residue_name));

		/* Translate and discretize the coordinates */
		cx = static_cast<uint16_t>((atm.x + ptran.x) * resolution + 0.5);
		cy = static_cast<uint16_t>((atm.y + ptran.y) * resolution + 0.5);
		cz = static_cast<uint16_t>((atm.z + ptran.z) * resolution + 0.5);

		radius = static_cast<uint16_t>((atm.radius + probeRadius) * resolution + 0.5);

		// x^2 + y^2 + z^2 = r^2
		for (int x = -radius; x <= radius; ++x) {
			// y^2 + z^2 = r^2 - x^2
			y_max = static_cast<int>(sqrt(radius * radius - x * x) + 0.5);
			for (int y = -y_max; y <= y_max; ++y) {
				// z^2 = r^2 - x^2 - y^2
				z_max = static_cast<int>(sqrt(radius * radius - x * x - y * y) + 0.5);
				for (int z = -z_max; z <= z_max; ++z) {
					if (surface->getVoxel(cx + x, cy + y, cz + z) && !(x == 0 && y == 0 && z == 0)) {
						float r = sqrt(x * x + y * y + z * z);

						BLAM930101_[cx + x][cy + y][cz + z] += hqi8[0] / r;
						BIOV880101_[cx + x][cy + y][cz + z] += hqi8[1] / r;
						MAXF760101_[cx + x][cy + y][cz + z] += hqi8[2] / r;
						TSAJ990101_[cx + x][cy + y][cz + z] += hqi8[3] / r;
						NAKH920108_[cx + x][cy + y][cz + z] += hqi8[4] / r;
						CEDJ970104_[cx + x][cy + y][cz + z] += hqi8[5] / r;
						LIFS790101_[cx + x][cy + y][cz + z] += hqi8[6] / r;
						MIYS990104_[cx + x][cy + y][cz + z] += hqi8[7] / r;
					}
				}
			}
		}
	}
}

