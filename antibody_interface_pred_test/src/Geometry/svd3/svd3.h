/*
 * svd3.h
 *
 *  Created on: 30/lug/2014
 *      Author: sebastian
 */

#ifndef DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_SVD3_SVD3_H_
#define DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_SVD3_SVD3_H_

#include "../eig3/eig3.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define ERR 1.0e-7
#define flt_equals(a, b) (fabs((a)-(b)) < ERR)


/**
 * Transposes the input 3x3 matrix and puts the resul in AT.
 */
static inline void transpose(double const A[3][3], double AT[3][3]) {
	AT[0][0] = A[0][0];
	AT[0][1] = A[1][0];
	AT[0][2] = A[2][0];
	AT[1][0] = A[0][1];
	AT[1][1] = A[1][1];
	AT[1][2] = A[2][1];
	AT[2][0] = A[0][2];
	AT[2][1] = A[1][2];
	AT[2][2] = A[2][2];
};
/**
 * Multiplies matrix A with matrix B and returns the result in C.
 */
static inline void matrix_multiply(double const A[3][3], double const B[3][3], double C[3][3]) {
	C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
	C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
	C[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];
	C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
	C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
	C[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];
	C[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
	C[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
	C[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
};

/**
 * This method calculates the Singular Value Decomposition of a 3x3 matrix.
 * A = S U V^T
 * @param A		matrix to decompose
 * @return U	orthogonal matrix, whose columns are eigenvectors of A A^T
 * @return S	diagonal matrix containing the singular values of A in decreasing order
 * @return V	orthogonal matrix, whose columns are eigenvectors of A^T A
 */
static inline void singular_value_decomposition(double const A[3][3],
		double U[3][3], double S[3][3], double V[3][3]) {
	/*
	 * Steps:
	 * 1) Use eigendecomposition on A^T A to compute V.
	 * Since A = U S V^T then A^T A = V S^T S V^T with D1 = S^T S and V the
	 * eigenvalues and eigenvectors respectively (V is orthogonal).
	 * 2) Compute U from A A^T, since A A^T=U S V^T V S^T U^T =U D2 U^T
	 * 3) Compute S as a diagonal 3x3 matrix with
	 */
	double AT[3][3], ATA[3][3], AAT[3][3], s1[3], s2[3];

	transpose(A, AT);

	matrix_multiply(AT, A, ATA);
	matrix_multiply(A, AT, AAT);

	eigen_decomposition(ATA, V, s1);
	eigen_decomposition(AAT, U, s2);

	assert(flt_equals(s1[0], s2[0])
			&& flt_equals(s1[1], s2[1])
			&& flt_equals(s1[2], s2[2]));

	S[0][0] = sqrt(s1[0]);
	S[0][1] = 0;
	S[0][2] = 0;
	S[1][0] = 0;
	S[1][1] = sqrt(s1[1]);
	S[1][2] = 0;
	S[2][0] = 0;
	S[2][1] = 0;
	S[2][2] = sqrt(s1[2]);
};
#endif /* DESCRIPTOREVALUATION_BACKUP_SRC_GEOMETRY_SVD3_SVD3_H_ */
