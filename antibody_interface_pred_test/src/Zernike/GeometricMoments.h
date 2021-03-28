
#ifndef GEOMETRICMOMENTS_H_
#define GEOMETRICMOMENTS_H_
/**
 @article{Hosny20071214,
 title = "Exact and fast computation of geometric moments for gray level images ",
 journal = "Applied Mathematics and Computation ",
 volume = "189",
 number = "2",
 pages = "1214 - 1222",
 year = "2007",
 issn = "0096-3003",
 doi = "http://dx.doi.org/10.1016/j.amc.2006.12.025",
 author = "Khalid M. Hosny"
 }
 */
#include "../Geometry/point3D.h"
#include <iostream>
#include <vector>

using namespace std;

typedef vector<double>				array1D;		// 1D array of scalar type
typedef vector<array1D>				array2D;		// 2D array of scalar type
typedef vector<array2D> 			array3D;        // 3D array of scalar type


typedef class GeometricMoments {
public:
#ifdef STATIC_GM
    static int maxOrder;       // maximal order of the moments
#else
    int maxOrder;
#endif
	/**
	 * Constructor of the GeometricMoments class.
	 * @param voxels	the input voxel grid
	 * @param scale		scale factor, in order to scale the given
	 * 					3D function inside the unit ball
	 * @param COG 		center of gravity of the input 3D function
	 * @param maxOrder		maximal order to compute moments for
	 */
	GeometricMoments(array3D const & voxels, double scale, point3D const & COG, int maxOrder);

	/**
	 * Null moments constructor
	 */
	GeometricMoments(int maxOrder);
	/**
	 * Default constructor.
	 */
	GeometricMoments();

	/**
	 * Returns the geometrical moment M_{rst} of order r+s+t,
	 * where r, s and t are the powers of x, y and z in the
	 * monomial x^r * y^s * z^t.
	 *
	 * @param r		order along x
	 * @param s		order along y
	 * @param t		order along z
	 * @return 		the value of the geometrical moment
	 */
	inline double getMoment(int r, int s, int t) const {
		return moments[r][s][t];
	};
private:
#ifdef STATIC_GM
	static int length, width, height;	// dimensions of the voxel grid
	static double dv; 	//size of the single voxel
	static point3D COG; 		// center of gravity
	static array1D x_samples, y_samples, z_samples; // scaled and translated samples
	static array2D Ir, Is, It;	// temporary data structures for the fast calculation
    							// of the geometric moments
#else
	int length, width, height;	// dimensions of the voxel grid
	double dv; 	//size of the single voxel
	point3D COG; 		// center of gravity
	array1D x_samples, y_samples, z_samples; // scaled and translated samples
	array2D Ir, Is, It;	// temporary data structures for the fast calculation
    					// of the geometric moments
#endif


    array3D const * voxels;	// pointer to the array containing the voxel grid

    array3D	moments;	// array containing the cumulative moments
	// ---- private methods ----
	void computeSamples();
	void computeI();
	void computeI_alt();
	void computeI_alt2();
	void computeI_old();

	void compute();
	// ---- debugging ----
	void compute_failsafe();
	void compute_failsafe_old();


} GeometricMoments;

#endif /* GEOMETRICMOMENTS_H_ */
