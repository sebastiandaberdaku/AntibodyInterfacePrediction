/*
 * point3D.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_POINT3D_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_POINT3D_H_

#include "../utils/doubleCompare.h"
#include <limits>
#include <math.h>
#include <stdexcept>
#include <ostream>

struct atom;

/**
 * Struct defining a 3D point data type.
 * A point3D identifies a single point in the 3D Cartesian space.
 */

typedef struct point3D {
private:
	static inline const double safeAcos(double const & x) {
		if (x < -1.0)
			return acos(-1.0);
		else if (x > 1.0)
			return acos(1.0);
		return acos(x);
	};
public:
	double x /**< x-coordinate */;
	double y /**< y-coordinate */;
	double z /**< z-coordinate */;

	//----------------
	//--Constructors--
	//----------------

	/**
	 * Default constructor.
	 * Initializes a point3D struct as (0, 0, 0).
	 */
	point3D() :
			x(0), y(0), z(0) { };
	/**
	 * Initializes a point3D struct as (x, y, z).
	 */
	point3D(double x, double y, double z) :
			x(x), y(y), z(z) { };
	/**
	 * Copy constructor.
	 */
	point3D(point3D const & p) :
			x(p.x), y(p.y), z(p.z) { };
	/**
	 * Copy assignment operator
	 */
	inline point3D & operator=(point3D const & p) {
		if (this != &p) {
			this->x = p.x;
			this->y = p.y;
			this->z = p.z;
		}
		return *this;
	};
	/**
	 * Constructs point as atom center.
	 */
	point3D(atom const & atm);
	//----------------------------------------------------------------
	//-- We use overloaded operators to implement vector operations.--
	//----------------------------------------------------------------

	/**
	 * Vector addition
	 */
	inline point3D & operator+=(point3D const & rhs) {
		this->x += rhs.x;
		this->y += rhs.y;
		this->z += rhs.z;
		return *this;
	};
	inline const point3D operator+(point3D const & rhs) const {
		return point3D(*this) += rhs;
	};
	/**
	 * Vector difference
	 */
	inline point3D & operator-=(point3D const & rhs) {
		this->x -= rhs.x;
		this->y -= rhs.y;
		this->z -= rhs.z;
		return *this;
	};
	inline const point3D operator-(point3D const & rhs) const {
		return point3D(*this) -= rhs;
	};
	inline const point3D operator-() const {
		return point3D(-(this->x), -(this->y), -(this->z));
	};
	/**
	 * Cross product
	 */
	inline point3D & operator*=(point3D const & rhs) {
		double x = this->y * rhs.z - this->z * rhs.y;
		double y = this->z * rhs.x - this->x * rhs.z;
		double z = this->x * rhs.y - this->y * rhs.x;
		this->x = x;
		this->y = y;
		this->z = z;
		return *this;
	};
	inline const point3D operator*(point3D const & rhs) const {
		return point3D(*this) *= rhs;
	};
	/**
	 * Product with scalar
	 */
	inline point3D & operator*=(double const & rhs) {
		this->x *= rhs;
		this->y *= rhs;
		this->z *= rhs;
		return *this;
	};
	inline const point3D operator*(double const & rhs) const {
		return point3D(*this) *= rhs;
	};
	friend inline const point3D operator*(double const & lhs, point3D const & rhs) {
		return point3D(rhs) *= lhs;
	};
	/**
	 * Transforms the current vector according to the transformation defined
	 * by the 3x3 matrix passed as argument.
	 */
	friend inline const point3D operator*(double lhs[3][3], point3D const & rhs) {
		double x = lhs[0][0] * rhs.x + lhs[0][1] * rhs.y + lhs[0][2] * rhs.z;
		double y = lhs[1][0] * rhs.x + lhs[1][1] * rhs.y + lhs[1][2] * rhs.z;
		double z = lhs[2][0] * rhs.x + lhs[2][1] * rhs.y + lhs[2][2] * rhs.z;
		return point3D(x, y, z);
	};
	/**
	 * Division with scalar
	 */
	inline point3D & operator/=(double const & rhs) {
		if (rhs == 0)
			throw std::invalid_argument("point3D::operator/=() - Divide by Zero!");
		this->x /= rhs;
		this->y /= rhs;
		this->z /= rhs;
		return *this;
	};
	inline const point3D operator/(double const & rhs) const {
		return point3D(*this) /= rhs;
	};
	/**
	 * Equal?
	 */
	inline bool operator==(point3D const & rhs) const {
		return doubleCompare(this->x, rhs.x)
				&& doubleCompare(this->y, rhs.y)
				&& doubleCompare(this->z, rhs.z);
	};
	/**
	 * Not equal?
	 */
	inline bool operator!=(point3D const & rhs) const {
		return !(*this == rhs);
	};
	/**
	 * Output to stream
	 */
	friend inline std::ostream & operator<<(std::ostream& os, point3D const & p) {
		os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
		return os;
	};

	//-----------
	//--Methods--
	//-----------
	/**
	 * Calculates the square of the distance between this point3D and the argument.
	 */
	inline double s_distance(point3D const & p) const {
		double dx = this->x - p.x;
		double dy = this->y - p.y;
		double dz = this->z - p.z;
		return (dx * dx + dy * dy + dz * dz);
	};
	/**
	 * Calculates the distance between this point3D and the argument.
	 */
	inline double distance(point3D const & p) const {
		return sqrt(this->s_distance(p));
	};
	/**
	 * Calculates the square of the distance between this point3D and the center of
	 * the argument.
	 */
	float s_distance(atom const & atm) const;
	/**
	 * Calculates the distance between this point3D and the center of the argument.
	 */
	inline float distance(atom const & atm) const {
		return sqrt(this->s_distance(atm));
	};
	/**
	 * Calculates the dot product between this point3D and the argument.
	 */
	inline const double dot(point3D const & p) const {
		return ((this->x) * p.x) + ((this->y) * p.y) + ((this->z) * p.z);
	};
	/**
	 * Calculates the square norm of the vector.
	 */
	inline const double s_norm() const{
		return this->dot(*this);
	};
	/**
	 * Calculates the norm of the vector.
	 */
	inline const double norm() const{
		return sqrt(this->s_norm());
	};

	/**
	 * Normalizes the current vector.
	 * If the current vector is the zero vector it does nothing.
	 */
	inline void normalize() {
		if (doubleCompare(this->x, 0) && doubleCompare(this->y, 0) && doubleCompare(this->z, 0))
			return;
		(*this) /= this->norm();
	};
	/**
	 * Calculates the angle (in radiant) between this point3D and the argument.
	 */
	inline const double angle(point3D const & p) const {
		return safeAcos((*this).dot(p)/(this->norm() * p.norm()));
	};
	/**
	 * Transforms the current vector according to the transformation defined
	 * by the 3x3 matrix passed as argument.
	 */
	inline const point3D transform(double const A[3][3]) const{
		double x = A[0][0] * this->x + A[0][1] * this->y + A[0][2] * this->z;
		double y = A[1][0] * this->x + A[1][1] * this->y + A[1][2] * this->z;
		double z = A[2][0] * this->x + A[2][1] * this->y + A[2][2] * this->z;
		return point3D(x, y, z);
	};
	/**
	 * Rotates the current vector by an angle of phi about an axis in the direction of u.
	 * @param u		Verse of the rotation (u must be a unitary vector)
	 * @param phi	Angle in radiant
	 *
	 * @return	The resulting vector.
	 */
	inline const point3D rotate(point3D const & u, double const & phi) const{
		if (!doubleCompare(1.0, u.norm()))
			throw std::invalid_argument("point3D::rotate() - Parameter u is not a unitary vector!");
		double R[3][3];
		R[0][0] = cos(phi) + u.x * u.x * (1 - cos(phi));
		R[0][1] = u.x * u.y * (1 - cos(phi)) - u.z * sin(phi);
		R[0][2] = u.x * u.z * (1 - cos(phi)) + u.y * sin(phi);
		R[1][0] = u.y * u.x * (1 - cos(phi)) + u.z * sin(phi);
		R[1][1] = cos(phi) + u.y * u.y * (1 - cos(phi));
		R[1][2] = u.y * u.z * (1 - cos(phi)) - u.x * sin(phi);
		R[2][0] = u.z * u.x * (1 - cos(phi)) - u.y * sin(phi);
		R[2][1] = u.z * u.y * (1 - cos(phi)) + u.x * sin(phi);
		R[2][2] = cos(phi) + u.z * u.z * (1 - cos(phi));
		return this->transform(R);
	};
} point3D;


#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_POINT3D_H_ */
