#pragma once
#include "corona_common.h" // for fast log and aligned alloc
#include "screenshot.h" // for fast log and aligned alloc
#include "filter.h"

#include <assert.h>
#include <string.h>
#include <math.h>

static inline void getEVDSymmetric2x2(double pOutEigenvalues[2], double pOutEigenvectors[4], const double pMatrix[4]) {
	// Define some short hands for the matrix entries
	double a = pMatrix[0];
	double b = pMatrix[1];
	double c = pMatrix[3];
	// Compute coefficients of the characteristic polynomial
	double pHalf = -0.5 * (a + c);
	double q = a*c - b*b;
	// Solve the quadratic
	double discriminant_root = sqrt(pHalf * pHalf - q);
	pOutEigenvalues[0] = -pHalf + discriminant_root;
	pOutEigenvalues[1] = -pHalf - discriminant_root;
	// Subtract a scaled identity matrix to obtain a rank one matrix
	double a0 = a - pOutEigenvalues[0];
	double b0 = b;
	double c0 = c - pOutEigenvalues[0];
	// The column space of this matrix is orthogonal to the first eigenvector 
	// and since the eigenvectors are orthogonal themselves, it agrees with the 
	// second eigenvector. Pick the longer column to avoid cancellation.
	double squaredLength0 = a0*a0 + b0*b0;
	double squaredLength1 = b0*b0 + c0*c0;
	double squaredLength;
	if (squaredLength0 > squaredLength1) {
		pOutEigenvectors[2] = a0;
		pOutEigenvectors[3] = b0;
		squaredLength = squaredLength0;
	}
	else {
		pOutEigenvectors[2] = b0;
		pOutEigenvectors[3] = c0;
		squaredLength = squaredLength1;
	}
	// If the eigenvector is exactly zero, both eigenvalues are the same and the 
	// choice of orthogonal eigenvectors is arbitrary
	pOutEigenvectors[2] = (squaredLength == 0.0) ? 1.0 : pOutEigenvectors[2];
	squaredLength = (squaredLength == 0.0) ? 1.0 : squaredLength;
	// Now normalize
	double invLength = 1.0 / sqrt(squaredLength);
	pOutEigenvectors[2] *= invLength;
	pOutEigenvectors[3] *= invLength;
	// And rotate to get the other eigenvector
	pOutEigenvectors[0] =  pOutEigenvectors[3];
	pOutEigenvectors[1] = -pOutEigenvectors[2];
}


/*! This utility function swaps the two doubles at the given pointers if they 
	are not sorted in descending order.*/
static inline void sortDescending2(double* pFirst, double* pSecond) {
	double temporary = *pSecond;
	int swap = *pFirst < temporary;
	(*pSecond) = swap ? *pFirst : temporary;
	(*pFirst) = swap ? temporary : (*pFirst);
}


/*! Sorts an array of three entries in descending order (in place).*/
static inline void sortDescending3(double pArray[3]) {
	sortDescending2(pArray + 0, pArray + 2);
	sortDescending2(pArray + 0, pArray + 1);
	sortDescending2(pArray + 1, pArray + 2);
}


/*! Given coefficients of a normalized, cubic polynomial
	x^3 + pCoefficient[2]*x^2  + pCoefficient[1]*x  + pCoefficient[0]
	that is known to have three real roots, this function computes and returns 
	all roots sorted in descending order. Behavior is undefined if there are 
	not three real roots.*/
static inline void solveCubic(double pOutRoots[3], const double pCoefficient[3]){
	// Divide middle coefficients by three
	double x = pCoefficient[0];
	double y = pCoefficient[1] / 3.0;
	double z = pCoefficient[2] / 3.0;
	// Compute the Hessian and the discrimant
	double pDelta[3] = {-z*z + y,-y*z + x,z*x - y*y};
	double discriminant = 4.0 * pDelta[0] * pDelta[2] - pDelta[1] * pDelta[1];
	double discriminantRoot = sqrt((discriminant < 0.0) ? 0.0 : discriminant);
	// Compute the root of largest magnitude
	double largestRoot;
	{
		// Compute coefficients of the depressed cubic (third is zero, fourth 
		// is one)
		double p = pDelta[0];
		double q = -2.0 * z * pDelta[0] + pDelta[1];
		// Take the cubic root of a normalized complex number (or maybe its 
		// conjugate, we only care about the real part)
		double theta = fabs(atan2(discriminantRoot, -q) / 3.0);
		double cubicRootX = cos(theta);
		double cubicRootY = sin(theta);
		// Take the real part of two of the three cubic roots. That literal is 
		// sqrt(3)/2.
		double rootCandidate0 = cubicRootX;
		double rootCandidate1 = -0.5 * cubicRootX - 0.86602540378443864676372317075294 * cubicRootY;
		// Apply the scaling of the depression transform
		double sqrtP = sqrt((p <= 0.0) ? (-p) : 0.0);
		rootCandidate0 *= 2.0 * sqrtP;
		rootCandidate1 *= 2.0 * sqrtP;
		// Pick the appropriate result
		largestRoot = ((rootCandidate0 + rootCandidate1) > 2.0 * z) ? rootCandidate0 : rootCandidate1;
		// Apply the shift of the depression transform
		largestRoot -= z;
	}
	// Compute the root of least magnitude
	double smallestRoot;
	{
		// Compute coefficients of the depressed cubic when working with a 
		// flipped polynomial where x is replaced by 1/x (third is zero)
		double A = x;
		double p = pDelta[2];
		double q = -A * pDelta[1] + 2.0 * y * pDelta[2];
		// Take the cubic root of a normalized complex number (or maybe its 
		// conjugate, we only care about the real part)
		double theta = fabs(atan2(A * discriminantRoot, -q) / 3.0);
		double cubicRootX = cos(theta);
		double cubicRootY = sin(theta);
		// Take the real part of two of the three cubic roots. That literal is 
		// sqrt(3)/2.
		double rootCandidate0 = cubicRootX;
		double rootCandidate1 = -0.5 * cubicRootX - 0.86602540378443864676372317075294 * cubicRootY;
		// Apply the scaling of the depression transform
		double sqrtP = sqrt((p <= 0.0) ? (-p) : 0.0);
		rootCandidate0 *= 2.0*sqrtP;
		rootCandidate1 *= 2.0*sqrtP;
		// Pick the appropriate result
		smallestRoot = (rootCandidate0 + rootCandidate1 < 2.0 * y) ? rootCandidate0 : rootCandidate1;
		// Apply the shift of the depression transform, account for the leading 
		// coefficient and undo the reciprocal transform
		smallestRoot = -A / (smallestRoot + y);
	}
	// Compute the root of intermediate magnitude through polynomial division
	double F = -largestRoot - smallestRoot;
	double G =  largestRoot * smallestRoot;
	double intermediateRoot = (y*F - z*G) / (-z*F + y);
	// Output and sort the roots (thus far they are only sorted by magnitude)
	pOutRoots[0] = largestRoot;
	pOutRoots[1] = intermediateRoot;
	pOutRoots[2] = smallestRoot;
	sortDescending3(pOutRoots);
}


/*! This utility function computes the cross-product of a vector and the 
	squared length of the cross product.*/
static inline void cross(double pOutCrossProduct[3], double lhs0, double lhs1, double lhs2, double rhs0, double rhs1, double rhs2) {
	pOutCrossProduct[0] = lhs1 * rhs2 - lhs2 * rhs1;
	pOutCrossProduct[1] = lhs2 * rhs0 - lhs0 * rhs2;
	pOutCrossProduct[2] = lhs0 * rhs1 - lhs1 * rhs0;
}

/*! Given a 3x3 matrix (or three 3D vectors, one after the other) this function 
	outputs a normalized version the column of greatest norm (or the longest 
	vector). If the matrix is zero entirely, it outputs an arbitrary normalized 
	vector.*/
static inline void selectLongestColumnAndNormalize(double pOutColumn[3], const double pMatrix[9]) {
	double maxSquaredLength = -1.0;
	int iLongest = 0;
	for (int i = 0; i != 3; ++i) {
		double squaredLength = pMatrix[3 * i + 0] * pMatrix[3 * i + 0] + pMatrix[3 * i + 1] * pMatrix[3 * i + 1] + pMatrix[3 * i + 2] * pMatrix[3 * i + 2];
		iLongest = (maxSquaredLength < squaredLength) ? i : iLongest;
		maxSquaredLength = (maxSquaredLength < squaredLength) ? squaredLength : maxSquaredLength;
	}
	if (maxSquaredLength == 0.0) {
		pOutColumn[0] = 1.0;
		pOutColumn[1] = 0.0;
		pOutColumn[2] = 0.0;
	}
	else {
		double invLength = 1.0 / sqrt(maxSquaredLength);
		for (int i = 0; i != 3; ++i) {
			pOutColumn[i] = invLength * pMatrix[3 * iLongest + i];
		}
	}
}


/*! This utility function checks whether the given two doubleing point values 
	are nearly equal.*/
static inline int equalWithTolerance(double lhs, double rhs) {
	const double epsilon = 1.0e-5f;
	int equalSign = ((lhs * rhs) >= 0.0);
	lhs = fabs(lhs);
	rhs = fabs(rhs);
	return equalSign && lhs <= (1.0 + epsilon) * rhs && lhs >= (1.0 - epsilon) * rhs;
}


static inline void getEVDSymmetric3x3(double pOutEigenvalues[3], double pOutEigenvectors[9], const double pMatrix[9]) {
	// Define some short hands for the matrix entries
	double a00 = pMatrix[0];
	double a10 = pMatrix[1];
	double a20 = pMatrix[2];
	double a11 = pMatrix[4];
	double a21 = pMatrix[5];
	double a22 = pMatrix[8];
	// Compute coefficients of the characteristic polynomial.
	double pCoefficient[3];
	// The constant one is minus the determinant
	pCoefficient[0] = -a00*(a11*a22 - a21 * a21) + a10*(a10*a22 - a20*a21) - a20*(a10*a21 - a11*a20);
	pCoefficient[1] = (a00*a11 + a11*a22 + a00*a22) - (a21*a21 + a20*a20 + a10*a10);
	// The quadratic one is minus the trace
	pCoefficient[2] = -a00 - a11 - a22;
	// Compute the eigenvalues
	solveCubic(pOutEigenvalues, pCoefficient);
	// The case of two nearly identical eigenvalues needs to be handled 
	// separately for robustness reasons
	int iUniqueEigenvalue = -1;
	int iSecondEigenvalue, iThirdEigenvalue;
	if (equalWithTolerance(pOutEigenvalues[0], pOutEigenvalues[1])) {
		iUniqueEigenvalue = 2;
		iSecondEigenvalue = 0;
		iThirdEigenvalue = 1;
	}
	else if (equalWithTolerance(pOutEigenvalues[1], pOutEigenvalues[2])) {
		iUniqueEigenvalue = 0;
		iSecondEigenvalue = 1;
		iThirdEigenvalue = 2;
	}
	if (iUniqueEigenvalue >= 0) {
		// Make the two eigenvalues truly equal
		pOutEigenvalues[iSecondEigenvalue] = pOutEigenvalues[iThirdEigenvalue] = 0.5 * (pOutEigenvalues[iSecondEigenvalue] + pOutEigenvalues[iThirdEigenvalue]);
		double doubleEigenvalue = pOutEigenvalues[2 - iUniqueEigenvalue];
		// The eigenvector for the unique eigenvalue is in the one-dimensional 
		// column space of the matrix minus the other eigenvalue
		double pSubtractedMatrix[9] = {
			a00 - doubleEigenvalue, a10                   , a20,
			a10                   , a11 - doubleEigenvalue, a21,
			a20                   , a21                   , a22 - doubleEigenvalue,
		};
		double pUniqueEigenvector[3];
		selectLongestColumnAndNormalize(pUniqueEigenvector, pSubtractedMatrix);
		// The other two vectors are orthonormal but arbitrary
		double pOrthogonal[9];
		cross(pOrthogonal + 0, pUniqueEigenvector[0], pUniqueEigenvector[1], pUniqueEigenvector[2], 1.0, 0.0, 0.0);
		cross(pOrthogonal + 3, pUniqueEigenvector[0], pUniqueEigenvector[1], pUniqueEigenvector[2], 0.0, 1.0, 0.0);
		cross(pOrthogonal + 6, pUniqueEigenvector[0], pUniqueEigenvector[1], pUniqueEigenvector[2], 0.0, 0.0, 1.0);
		double pSecondEigenvector[3];
		selectLongestColumnAndNormalize(pSecondEigenvector, pOrthogonal);
		double pThirdEigenvector[3];
		cross(pThirdEigenvector, pUniqueEigenvector[0], pUniqueEigenvector[1], pUniqueEigenvector[2], pSecondEigenvector[0], pSecondEigenvector[1], pSecondEigenvector[2]);
		// Store the eigenvectors
		for (int i = 0; i != 3; ++i) {
			pOutEigenvectors[3 * iUniqueEigenvalue + i] = pUniqueEigenvector[i];
			pOutEigenvectors[3 * iSecondEigenvalue + i] = pSecondEigenvector[i];
			pOutEigenvectors[3 * iThirdEigenvalue  + i] = pThirdEigenvector[i];
		}
		return;
	}
	// Otherwise, we compute the first two eigenvectors independently
	for (int i = 0; i != 2; ++i) {
		// Subtract the identity matrix times the eigenvalue
		double b00 = a00 - pOutEigenvalues[i];
		double b10 = a10;
		double b20 = a20;
		double b11 = a11 - pOutEigenvalues[i];
		double b21 = a21;
		double b22 = a22 - pOutEigenvalues[i];
		// The eigenvector is orthogonal to the column space of this matrix. 
		// Find it by taking the cross product of two columns. For better 
		// stabiliy, we try all three pairs.
		double pCandidate[9];
		cross(pCandidate + 0, b00, b10, b20, b10, b11, b21);
		cross(pCandidate + 3, b10, b11, b21, b20, b21, b22);
		cross(pCandidate + 6, b20, b21, b22, b00, b10, b20);
		selectLongestColumnAndNormalize(pOutEigenvectors + 3 * i, pCandidate);
		// Explicitly enforce orthogonality
		if (i == 1) {
			double dot = pOutEigenvectors[0] * pOutEigenvectors[3] + pOutEigenvectors[1] * pOutEigenvectors[4] + pOutEigenvectors[2] * pOutEigenvectors[5];
			pOutEigenvectors[3] -= dot * pOutEigenvectors[0];
			pOutEigenvectors[4] -= dot * pOutEigenvectors[1];
			pOutEigenvectors[5] -= dot * pOutEigenvectors[2];
			// Renormalize
			double invLength = 1.0 / sqrt(pOutEigenvectors[3] * pOutEigenvectors[3] + pOutEigenvectors[4] * pOutEigenvectors[4] + pOutEigenvectors[5] * pOutEigenvectors[5]);
			pOutEigenvectors[3] *= invLength;
			pOutEigenvectors[4] *= invLength;
			pOutEigenvectors[5] *= invLength;
		}
	}
	// The last one is orthogonal to the first two
	cross(pOutEigenvectors + 6, pOutEigenvectors[0], pOutEigenvectors[1], pOutEigenvectors[2], pOutEigenvectors[3], pOutEigenvectors[4], pOutEigenvectors[5]);
}
