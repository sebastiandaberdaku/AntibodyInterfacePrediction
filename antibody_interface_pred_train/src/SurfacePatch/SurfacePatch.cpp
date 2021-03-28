/*
 * SurfacePatch.cpp
 *
 *  Created on: Jul 16, 2014
 *      Author: sebastian
 */
/**
 * Implementation of the SurfacePatch object.
 */
#include "../Geometry/normalizeGrid.h"
#include "../Geometry/Sphere/SphereOffsets.h"
#include "../utils/makeDirectory.h"
#include "../Zernike/BoundingSphere.h"
#include "SurfacePatch.h"
#include <iomanip>
#include <queue>

using namespace std;

list<voxel_offset> SurfacePatch::patchSphereOffsets;

SurfacePatch::SurfacePatch(size_t id, float patchRadius,
		voxel const & center, float resolution,
		voxelGrid const * molecularSurface, voxelGrid const * volumetricModel,
		array3D const & potentials,
		array3D const & hydrophobicity,
		point3D const & ptran, array3D const & interface,
		array3D const & _BLAM930101,
		array3D const & _BIOV880101,
		array3D const & _MAXF760101,
		array3D const & _TSAJ990101,
		array3D const & _NAKH920108,
		array3D const & _CEDJ970104,
		array3D const & _LIFS790101,
		array3D const & _MIYS990104) :
		ID(id), patchRadius(patchRadius), patchCenter(center), surfaceD(NULL),
		potentials_posD(NULL), potentials_negD(NULL),
		hydrophobicity_posD(NULL), hydrophobicity_negD(NULL),
		BLAM930101_posD(NULL),
		BLAM930101_negD(NULL),
		BIOV880101_posD(NULL),
		BIOV880101_negD(NULL),
		MAXF760101_D(NULL),
		TSAJ990101_D(NULL),
		NAKH920108_D(NULL),
		CEDJ970104_D(NULL),
		LIFS790101_D(NULL),
		MIYS990104_posD(NULL),
		MIYS990104_negD(NULL),
		translation(ptran), resolution(resolution),
		molecularSurface(molecularSurface), volumetricModel(volumetricModel) {

	int x, y, z;
	int d_r = static_cast<int>(patchRadius * resolution + 0.5);
	dim = 2 * d_r + 1;
#pragma omp critical
{
	if (patchSphereOffsets.empty())
		get_sphere_offsets(d_r, patchSphereOffsets);
}
//	this->surface            = array3D(dim, array2D(dim, array1D(dim, 0)));
//	this->hydrophobicity_pos = array3D(dim, array2D(dim, array1D(dim, 0)));
//	this->hydrophobicity_neg = array3D(dim, array2D(dim, array1D(dim, 0)));
//	this->potentials_pos = array3D(dim, array2D(dim, array1D(dim, 0)));
//	this->potentials_neg = array3D(dim, array2D(dim, array1D(dim, 0)));

	BLAM930101_pos = array3D(dim, array2D(dim, array1D(dim, 0)));
	BLAM930101_neg = array3D(dim, array2D(dim, array1D(dim, 0)));
	BIOV880101_pos = array3D(dim, array2D(dim, array1D(dim, 0)));
	BIOV880101_neg = array3D(dim, array2D(dim, array1D(dim, 0)));
	MAXF760101_    = array3D(dim, array2D(dim, array1D(dim, 0)));
	TSAJ990101_    = array3D(dim, array2D(dim, array1D(dim, 0)));
	NAKH920108_    = array3D(dim, array2D(dim, array1D(dim, 0)));
	CEDJ970104_    = array3D(dim, array2D(dim, array1D(dim, 0)));
	LIFS790101_    = array3D(dim, array2D(dim, array1D(dim, 0)));
	MIYS990104_pos = array3D(dim, array2D(dim, array1D(dim, 0)));
	MIYS990104_neg = array3D(dim, array2D(dim, array1D(dim, 0)));

	size_t patch_area = 0;
	size_t interface_area = 0;
	size_t patch_volume = 0;

	for (auto const & o : patchSphereOffsets) {
		x = center.ix + o.i;
		y = center.iy + o.j;
		z = center.iz + o.k;
		if (x >= 0 && x < molecularSurface->length
				&& y >= 0 && y < molecularSurface->width
				&& z >= 0 && z < molecularSurface->height){
			if (molecularSurface->getVoxel(x, y, z)) {
//				this->surface[o.i + d_r][o.j + d_r][o.k + d_r] = 1;
				++patch_area;
				++patch_volume;
//				if (potentials[x][y][z] > 0)
//					potentials_pos[o.i + d_r][o.j + d_r][o.k + d_r] = potentials[x][y][z];
//				else
//					potentials_neg[o.i + d_r][o.j + d_r][o.k + d_r] = -potentials[x][y][z];
//				if (hydrophobicity[x][y][z] > 0)
//					this->hydrophobicity_pos[o.i + d_r][o.j + d_r][o.k + d_r] = hydrophobicity[x][y][z];
//				else
//					this->hydrophobicity_neg[o.i + d_r][o.j + d_r][o.k + d_r] = -hydrophobicity[x][y][z];

				if (_BLAM930101[x][y][z] > 0)
					BLAM930101_pos[o.i + d_r][o.j + d_r][o.k + d_r] = _BLAM930101[x][y][z];
				else
					BLAM930101_neg[o.i + d_r][o.j + d_r][o.k + d_r] = -_BLAM930101[x][y][z];
				if (_BIOV880101[x][y][z] > 0)
					BIOV880101_pos[o.i + d_r][o.j + d_r][o.k + d_r] = _BIOV880101[x][y][z];
				else
					BIOV880101_neg[o.i + d_r][o.j + d_r][o.k + d_r] = -_BIOV880101[x][y][z];
				MAXF760101_[o.i + d_r][o.j + d_r][o.k + d_r]    = _MAXF760101[x][y][z];
				TSAJ990101_[o.i + d_r][o.j + d_r][o.k + d_r]    = _TSAJ990101[x][y][z];
				NAKH920108_[o.i + d_r][o.j + d_r][o.k + d_r]    = _NAKH920108[x][y][z];
				CEDJ970104_[o.i + d_r][o.j + d_r][o.k + d_r]    = _CEDJ970104[x][y][z];
				LIFS790101_[o.i + d_r][o.j + d_r][o.k + d_r]    = _LIFS790101[x][y][z];
				if (_MIYS990104[x][y][z] > 0)
					MIYS990104_pos[o.i + d_r][o.j + d_r][o.k + d_r] = _MIYS990104[x][y][z];
				else
					MIYS990104_neg[o.i + d_r][o.j + d_r][o.k + d_r] = -_MIYS990104[x][y][z];
			}
			else if (volumetricModel->getVoxel(x, y, z)) {
				++patch_volume;
			}
			if (interface[x][y][z] > 0) {
				++interface_area;
			}
		}
	}

	this->interface = (interface[center.ix][center.iy][center.iz] > 0);
//	this->interface = (interfaceThreshold <= (interface_area / (float) patch_area));
//	this->curvature = patch_volume / (float) patchSphereOffsets.size(); // is the sphere volume in voxels
//
//
//
//	int32_t nx, ny, nz;
//
//	double C[3][3]; /**> correlation matrix of the current patch */
//	double eigenvectors[3][3]; /**> each column will be an eigenvector of C*/
//	double eigenvalues[3]; /**> eigenvalues corresponding to the eigenvectors*/
//
//	calculateCorrelationMatrix(C);
//	eigen_decomposition(C, eigenvectors, eigenvalues);
//	/** The normal can be approximated with the eigenvector corresponding
//	 * to the smallest eigenvalue. */
//	patchNormal.x = eigenvectors[0][2];
//	patchNormal.y = eigenvectors[1][2];
//	patchNormal.z = eigenvectors[2][2];
//
//	/** Now we determine the normal's orientation */
//	nx = patchCenter.ix + (d_r / 2) * patchNormal.x;
//	ny = patchCenter.iy + (d_r / 2) * patchNormal.y;
//	nz = patchCenter.iz + (d_r / 2) * patchNormal.z;
//
//	if (nx >= 0 && nx < volumetricModel->length && ny >= 0
//			&& ny < volumetricModel->width && nz >= 0
//			&& nz < volumetricModel->height) {
//		if (volumetricModel->getVoxel(nx, ny, nz)) {
//			patchNormal *= -1;
//		}
//	}
//	/** Clearly if going out of bounds, the normal points
//	 * in the correct direction. */
//	/**
//	 * flatness is a real between 0 and 1/3
//	 */
//	flatness = eigenvalues[2]
//			/ (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
}

SurfacePatch::~SurfacePatch() {
//	if (surfaceD != NULL)
//		delete surfaceD;
//	if (potentials_posD != NULL)
//		delete potentials_posD;
//	if (potentials_negD != NULL)
//		delete potentials_negD;
//	if (hydrophobicity_posD != NULL)
//		delete hydrophobicity_posD;
//	if (hydrophobicity_negD != NULL)
//		delete hydrophobicity_negD;
	if (BLAM930101_posD != NULL)
		delete BLAM930101_posD;
	if (BLAM930101_negD != NULL)
		delete BLAM930101_negD;
	if (BIOV880101_posD != NULL)
		delete BIOV880101_posD;
	if (BIOV880101_negD != NULL)
		delete BIOV880101_negD;
	if (MAXF760101_D  != NULL)
		delete MAXF760101_D;
	if (TSAJ990101_D != NULL)
		delete TSAJ990101_D;
	if (NAKH920108_D != NULL)
		delete NAKH920108_D;
	if (CEDJ970104_D != NULL)
		delete CEDJ970104_D;
	if (LIFS790101_D != NULL)
		delete LIFS790101_D;
	if (MIYS990104_posD != NULL)
		delete MIYS990104_posD;
	if (MIYS990104_negD  != NULL)
		delete MIYS990104_negD;
}

void SurfacePatch::clearDisconnectedVoxels() {
	int c = dim / 2; // notice that dim is always odd

	assert(surface[c][c][c] > 0);

	int nbX, nbY, nbZ;

	array3D temp_s(surface);
	array3D temp_p(potentials_pos);
	array3D temp_n(potentials_neg);

	surface = array3D(dim, array2D(dim, array1D(dim, 0)));
	potentials_pos = array3D(dim, array2D(dim, array1D(dim, 0)));
	potentials_neg = array3D(dim, array2D(dim, array1D(dim, 0)));

	std::queue<voxel> pendingVoxels;
	voxel cv(c, c, c); // current voxel
	pendingVoxels.push(cv);
	temp_s[cv.ix][cv.iy][cv.iz] = 0;
	while (!pendingVoxels.empty()) {
		cv = pendingVoxels.front();
		pendingVoxels.pop();
		surface[cv.ix][cv.iy][cv.iz] = 1;
		potentials_pos[cv.ix][cv.iy][cv.iz] = temp_p[cv.ix][cv.iy][cv.iz];
		potentials_neg[cv.ix][cv.iy][cv.iz] = temp_n[cv.ix][cv.iy][cv.iz];
		for (int i = 0; i < 26; ++i) {
			/* try to fill currentVoxel's neighbors */
			nbX = cv.ix + voxel_offset::neighbours[i].i;
			nbY = cv.iy + voxel_offset::neighbours[i].j;
			nbZ = cv.iz + voxel_offset::neighbours[i].k;
			/* If not going out of bounds (i.e. if the neighbor exists),
			 * if the current neighbor has not been colored yet, and if
			 * it is not occupied put it in the queue. */
			if (nbZ > -1 && nbZ < dim && nbY > -1 && nbY < dim && nbX > -1 && nbX < dim) {
				if (temp_s[nbX][nbY][nbZ] > 0 && surface[nbX][nbY][nbZ] == 0) {
					// color it!
					temp_s[nbX][nbY][nbZ] = 0;
					pendingVoxels.push(voxel(nbX, nbY, nbZ));
				} // if
			} // if
		} // for
	} // while
}

void SurfacePatch::normalEstimator() {
	int32_t nx, ny, nz;
	/* discretized patch radius */

	double C[3][3]; /**> correlation matrix of the current patch */
	double eigenvectors[3][3]; /**> each column will be an eigenvector of C*/
	double eigenvalues[3]; /**> eigenvalues corresponding to the eigenvectors*/

	calculateCorrelationMatrix(C);
	eigen_decomposition(C, eigenvectors, eigenvalues);
	/** The normal can be approximated with the eigenvector corresponding
	 * to the smallest eigenvalue. */
	patchNormal.x = eigenvectors[0][2];
	patchNormal.y = eigenvectors[1][2];
	patchNormal.z = eigenvectors[2][2];

	int d_radius = dim / 2; // dim is odd

	/** Now we determine the normal's orientation */
	nx = patchCenter.ix + (d_radius / 2) * patchNormal.x;
	ny = patchCenter.iy + (d_radius / 2) * patchNormal.y;
	nz = patchCenter.iz + (d_radius / 2) * patchNormal.z;

	if (nx >= 0 && nx < volumetricModel->length && ny >= 0
			&& ny < volumetricModel->width && nz >= 0
			&& nz < volumetricModel->height) {
		if (volumetricModel->getVoxel(nx, ny, nz)) {
			patchNormal *= -1;
		}
	}
	/** Clearly if going out of bounds, the normal points
	 * in the correct direction. */
//	float curvature = eigenvalues[2]
//			/ (eigenvalues[0] + eigenvalues[1] + eigenvalues[2]);
//	/**
//	 * curvature is a real between 0 and 1/3
//	 */
//	if (curvature < alpha) {
//		type = 0;
//	} else {
//		int COGx = static_cast<int>(ceil(COG.x + 0.5)) + patchCenter.ix - d_radius;
//		int COGy = static_cast<int>(ceil(COG.y + 0.5)) + patchCenter.iy - d_radius;
//		int COGz = static_cast<int>(ceil(COG.z + 0.5)) + patchCenter.iz - d_radius;
//		if (volumetricModel->getVoxel(COGx, COGy, COGz))
//			type = 1;
//		else
//			type = -1;
//	}
}

void SurfacePatch::calculateZernikeDescriptors(int order) {
#ifdef STATIC_GM
	// dim is always odd
	int r = dim / 2;
	point3D center(r, r, r);
	double scale = 1.0 / (r + 1);
#else
	point3D center = calculateCOG(surface);
	double scale = 1.0 / (buondingSphereRadius(surface, center) + 1);
#endif
//	surfaceD = new ZernikeDescriptor(surface, center, scale, order);
//	potentials_posD = new ZernikeDescriptor(potentials_pos, center, scale, order);
//	potentials_negD = new ZernikeDescriptor(potentials_neg, center, scale, order);
//	hydrophobicity_posD = new ZernikeDescriptor(hydrophobicity_pos,     center, scale, order);
//	hydrophobicity_negD = new ZernikeDescriptor(hydrophobicity_neg,     center, scale, order);

	BLAM930101_posD = new ZernikeDescriptor(BLAM930101_pos, center, scale, order);
	BLAM930101_negD = new ZernikeDescriptor(BLAM930101_neg, center, scale, order);
	BIOV880101_posD = new ZernikeDescriptor(BIOV880101_pos, center, scale, order);
	BIOV880101_negD = new ZernikeDescriptor(BIOV880101_neg, center, scale, order);
	MAXF760101_D    = new ZernikeDescriptor(MAXF760101_, 	center, scale, order);
	TSAJ990101_D    = new ZernikeDescriptor(TSAJ990101_, 	center, scale, order);
	NAKH920108_D    = new ZernikeDescriptor(NAKH920108_, 	center, scale, order);
	CEDJ970104_D    = new ZernikeDescriptor(CEDJ970104_, 	center, scale, order);
	LIFS790101_D    = new ZernikeDescriptor(LIFS790101_, 	center, scale, order);
	MIYS990104_posD = new ZernikeDescriptor(MIYS990104_pos, center, scale, order);
	MIYS990104_negD = new ZernikeDescriptor(MIYS990104_neg, center, scale, order);
}

void SurfacePatch::exportPatchPotentialsPCDModel(string const & filename) {
	assert(!surface.empty());

	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open(("./output/" + filename + "-patch-" + to_string(ID) + ".pcd").c_str());

	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			for (int k = 0; k < dim; ++k) {
				if(surface[i][j][k] > 0)
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}
	union {
		uint8_t rgb_byte[4];
		float rgb_flt;
	} color_RGB;

	file_stream  << std::setprecision(11);
	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n" << "FIELDS x y z rgb\n" << "SIZE 4 4 4 "<< sizeof(color_RGB) <<"\n"
			<< "TYPE F F F F\n" << "COUNT 1 1 1 1\n" << "WIDTH "
			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
			<< "\n" << "DATA ascii";
	float pDX; /**< current potential */
	int x, y, z;
	uint8_t r, g, b;

	while (!surfaceVoxels.empty()) {
		x = surfaceVoxels.front().ix;
		y = surfaceVoxels.front().iy;
		z = surfaceVoxels.front().iz;

		pDX = potentials_pos[x][y][z] - potentials_neg[x][y][z];
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

		file_stream << "\n"
				<< surfaceVoxels.front().ix << " "
				<< surfaceVoxels.front().iy << " "
				<< surfaceVoxels.front().iz << " "
				<< color_RGB.rgb_flt;
		surfaceVoxels.pop();
	}
	file_stream.close();
}
void SurfacePatch::exportPatchPPotentialsPCDModel(string const & filename) {
	assert(!surface.empty());

	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open(("./output/" + filename + "-patch-" + to_string(ID) + ".pcd").c_str());

	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			for (int k = 0; k < dim; ++k) {
				if(surface[i][j][k] > 0)
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}
	union {
		uint8_t rgb_byte[4];
		float rgb_flt;
	} color_RGB;

	file_stream  << std::setprecision(11);
	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n" << "FIELDS x y z rgb\n" << "SIZE 4 4 4 "<< sizeof(color_RGB) <<"\n"
			<< "TYPE F F F F\n" << "COUNT 1 1 1 1\n" << "WIDTH "
			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
			<< "\n" << "DATA ascii";
	float pDX; /**< current potential */
	int x, y, z;
	uint8_t r, g, b;

	while (!surfaceVoxels.empty()) {
		x = surfaceVoxels.front().ix;
		y = surfaceVoxels.front().iy;
		z = surfaceVoxels.front().iz;

		pDX = potentials_pos[x][y][z];
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

		file_stream << "\n"
				<< surfaceVoxels.front().ix << " "
				<< surfaceVoxels.front().iy << " "
				<< surfaceVoxels.front().iz << " "
				<< color_RGB.rgb_flt;
		surfaceVoxels.pop();
	}
	file_stream.close();
}
void SurfacePatch::exportPatchNPotentialsPCDModel(string const & filename) {
	assert(!surface.empty());

	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open(("./output/" + filename + "-patch-" + to_string(ID) + ".pcd").c_str());

	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			for (int k = 0; k < dim; ++k) {
				if(surface[i][j][k] > 0)
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}
	union {
		uint8_t rgb_byte[4];
		float rgb_flt;
	} color_RGB;

	file_stream  << std::setprecision(11);
	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n" << "FIELDS x y z rgb\n" << "SIZE 4 4 4 "<< sizeof(color_RGB) <<"\n"
			<< "TYPE F F F F\n" << "COUNT 1 1 1 1\n" << "WIDTH "
			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
			<< "\n" << "DATA ascii";
	float pDX; /**< current potential */
	int x, y, z;
	uint8_t r, g, b;

	while (!surfaceVoxels.empty()) {
		x = surfaceVoxels.front().ix;
		y = surfaceVoxels.front().iy;
		z = surfaceVoxels.front().iz;

		pDX = - potentials_neg[x][y][z];
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

		file_stream << "\n"
				<< surfaceVoxels.front().ix << " "
				<< surfaceVoxels.front().iy << " "
				<< surfaceVoxels.front().iz << " "
				<< color_RGB.rgb_flt;
		surfaceVoxels.pop();
	}
	file_stream.close();
}

void SurfacePatch::exportPatchPCDModel(string const & filename) {
	assert(!surface.empty());

	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open(("./output/" + filename + "-patch-" + to_string(ID) + ".pcd").c_str());

	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			for (int k = 0; k < dim; ++k) {
				if(surface[i][j][k] > 0)
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}

	file_stream  << std::setprecision(11);

	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n"
			<< "FIELDS x y z\n"
			<< "SIZE 4 4 4\n" << "TYPE F F F\n"
			<< "COUNT 1 1 1\n" << "WIDTH " << surfaceVoxels.size()
			<< "\n" << "HEIGHT 1\n" << "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS "
			<< surfaceVoxels.size() << "\n" << "DATA ascii";
	int c = dim / 2; // notice that dim is always odd

	while (!surfaceVoxels.empty()) {
		file_stream << "\n" << (patchCenter.ix + surfaceVoxels.front().ix - c) /resolution - translation.x
				<< " " << (patchCenter.iy + surfaceVoxels.front().iy - c) /resolution - translation.y
				<< " " << (patchCenter.iz + surfaceVoxels.front().iz - c) /resolution - translation.z;
		surfaceVoxels.pop();
	}
	file_stream.close();
}

void SurfacePatch::exportZernikeInvariants(string const & filename) {
	assert(surfaceD != NULL);
	surfaceD->saveInvariants(filename + "-surface");
	if (potentials_posD != NULL)
		potentials_posD->saveInvariants(filename + "-pos_potentials");
	if (potentials_negD != NULL)
		potentials_negD->saveInvariants(filename + "-neg_potentials");
}

void SurfacePatch::reconstructPatch() {
	complexArray3D grid;
	    grid.resize(dim);
	    for (int i = 0; i < dim; ++i)
	    	grid[i].resize(dim);
	    for (int i = 0; i < dim; ++i)
	        for (int j = 0; j < dim; ++j)
	    	grid[i][j].resize(dim);

	    surfaceD->reconstruct(grid);
	    reconstructedPatch.resize(dim);

	for (int i = 0; i < dim; ++i) {
		reconstructedPatch[i].resize(dim);
		for (int j = 0; j < dim; ++j) {
			reconstructedPatch[i][j].resize(dim);
			for (int k = 0; k < dim; ++k) {
				reconstructedPatch[i][j][k] = norm(grid[i][j][k]);
			}
		}
	}
}

void SurfacePatch::exportReconstructedPatchPCDModel(string const & filename, double threshold) {
	assert(!reconstructedPatch.empty());

	std::ofstream file_stream;
	makeDirectory("./output");
	file_stream.open(("./output/" + filename + ".pcd").c_str());

	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < dim; ++j) {
			for (int k = 0; k < dim; ++k) {
				if(reconstructedPatch[i][j][k] > threshold)
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}
	file_stream  << std::setprecision(11);
	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "VERSION .7\n" << "FIELDS x y z\n" << "SIZE 4 4 4\n"
			<< "TYPE F F F\n" << "COUNT 1 1 1\n" << "WIDTH "
			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
			<< "\n" << "DATA ascii";

	while (!surfaceVoxels.empty()) {
		file_stream << "\n" << surfaceVoxels.front().ix
				<< " " << surfaceVoxels.front().iy
				<< " " << surfaceVoxels.front().iz;
		surfaceVoxels.pop();
	}
	file_stream.close();
}
