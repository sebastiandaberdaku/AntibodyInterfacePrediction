/*
 * PotentialGridDX.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_POTENTIALGRIDDX_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_POTENTIALGRIDDX_H_

#include "../Geometry/point3D.h"
#include <string>

using namespace std;

#define VGRID_DIGITS 6

typedef struct PotentialGridDX {
private:
    double Vcompare;
    /**
     * Returns the value at grid point [i, j, k].
     */
    inline double getValue(uint32_t i, uint32_t j, uint32_t k) {
    	return data[i * ny * nz + j * nz + k];
    };
public:
    int nx;       /**< Number of grid points in x direction */
    int ny;       /**< Number of grid points in y direction */
    int nz;       /**< Number of grid points in z direction */
    double hx;    /**< Grid spacing in x direction */
    double hy;    /**< Grid spacing in y direction */
    double hz;	  /**< Grid spacing in z direction */
    double xmin;  /**< x coordinate of lower grid corner */
    double ymin;  /**< y coordinate of lower grid corner */
    double zmin;  /**< z coordinate of lower grid corner */
    double xmax;  /**< x coordinate of upper grid corner */
    double ymax;  /**< y coordinate of upper grid corner */
    double zmax;  /**< z coordinate of upper grid corner */
    double *data; /**< nx*ny*nz array of data */
    PotentialGridDX(string const & openDX_Filename);
    virtual ~PotentialGridDX();
    /**
     * Returns the value of the electrostatic potential at a given 3D position.
     * If the requested point matches one of the grid points, then the corresponding
     * grid value is returned. Else, the value of the electrostatic potential is
     * calculated by trilinear interpolation.
     * If the requested 3D point falls outside the potentials grid, then a value
     * of 0 is returned.
     * @param p		position in space where we want to know the electrostatic
     * 				potential
     *
     * @return		the value of the electrostatic potential in the requested point
     * 				in space
     */
	inline double getPotential(point3D p) {
		int ihi, jhi, khi, ilo, jlo, klo;
		double ifloat, jfloat, kfloat;
		double u, dx, dy, dz;

		u = 0;

		ifloat = (p.x - xmin) / hx;
		jfloat = (p.y - ymin) / hy;
		kfloat = (p.z - zmin) / hz;

		ihi = (int) ceil(ifloat);
		jhi = (int) ceil(jfloat);
		khi = (int) ceil(kfloat);
		ilo = (int) floor(ifloat);
		jlo = (int) floor(jfloat);
		klo = (int) floor(kfloat);
		if (fabs(p.x - xmin) < Vcompare)
			ilo = 0;
		if (fabs(p.y - ymin) < Vcompare)
			jlo = 0;
		if (fabs(p.z - zmin) < Vcompare)
			klo = 0;
		if (fabs(p.x - xmax) < Vcompare)
			ihi = nx - 1;
		if (fabs(p.y - ymax) < Vcompare)
			jhi = ny - 1;
		if (fabs(p.z - zmax) < Vcompare)
			khi = nz - 1;

		/* See if we're on the mesh */
		if ((ihi < nx) && (jhi < ny) && (khi < nz) && (ilo >= 0) && (jlo >= 0)
				&& (klo >= 0)) {

			dx = ifloat - (double) (ilo);
			dy = jfloat - (double) (jlo);
			dz = kfloat - (double) (klo);
			u = dx * dy * dz * (getValue(ihi, jhi, khi))
					+ dx * (1.0 - dy) * dz * (getValue(ihi, jlo, khi))
					+ dx * dy * (1.0 - dz) * (getValue(ihi, jhi, klo))
					+ dx * (1.0 - dy) * (1.0 - dz) * (getValue(ihi, jlo, klo))
					+ (1.0 - dx) * dy * dz * (getValue(ilo, jhi, khi))
					+ (1.0 - dx) * (1.0 - dy) * dz * (getValue(ilo, jlo, khi))
					+ (1.0 - dx) * dy * (1.0 - dz) * (getValue(ilo, jhi, klo))
					+ (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
							* (getValue(ilo, jlo, klo));

			return u;

		} else
			return 0;
	};
} PotentialGridDX;

#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_POTENTIALGRIDDX_H_ */
