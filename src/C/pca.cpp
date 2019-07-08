// autotools
//#include <config.h>

// c++ header files
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <vector>

// gsl header files
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cdf.h>

// local includes
#include "common.hpp"

// namespaces
using namespace std;

//Throughout this module, a matrix of points has
//each row comprising one point and each column comprising one dimension.

typedef struct eigen_t {
	gsl_matrix *vectors;
	gsl_vector *values;
} Eigen;

void freeEigen(Eigen *eig) {
	gsl_matrix_free(eig->vectors);
	gsl_vector_free(eig->values);
};

//Find the eigenvalues and eigenvectors of a symmetric matrix, such as a covariance matrix.
//Returns them in descending order of eigenvalue.
//Eigenvectors are in columns in the matrix.
Eigen *getEigen(const gsl_matrix *mat) {
	if(mat->size1 != mat->size2) {
		fprintf(stderr, "Cannot calculate eigenvalues for a non-square matrix.\n");
		exit(1);
	}
	//Allocate the return struct.
	Eigen *eig = (Eigen*)malloc(sizeof(Eigen));
	eig->values = gsl_vector_alloc(mat->size1);
	eig->vectors = gsl_matrix_alloc(mat->size1, mat->size1);
	//The function needs additional workspace.
	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(mat->size1);
	//The function will destroy the matrix, which is const, so grab a copy of it.
	gsl_matrix *matCopy = gsl_matrix_alloc(mat->size1, mat->size2);
	gsl_matrix_memcpy(matCopy, mat);
	//Find and sort the eigenthings.
	gsl_eigen_symmv(matCopy, eig->values, eig->vectors, w);
	gsl_eigen_symmv_sort(eig->values, eig->vectors, GSL_EIGEN_SORT_VAL_DESC);
	
	//Free the workspace and the copy.
	gsl_eigen_symmv_free(w);
	gsl_matrix_free(matCopy);
	
	return eig;
}

//This function should only be necessary for debugging.
void printMat(const gsl_matrix *mat) {
	for(size_t i = 0; i < mat->size1; i++) {
		for(size_t j = 0; j < mat->size2; j++) {
			printf("%11.6lf", gsl_matrix_get(mat, i, j));
		}
		printf("\n");
	}
}

//Returns the Euclidean distance between two points
//in whichever co-ordinate system they are provided in.
double getDist(gsl_vector *a, gsl_vector *b) {
	if(a->size != b->size) {
		fprintf(stderr, "Cannot calculate the distance between two points in spaces of different dimensions.\n");
		exit(1);
	}
	double dist = 0.0;
	for(size_t i = 0; i < a->size; i++) {
		dist += pow(gsl_vector_get(a, i) - gsl_vector_get(b, i), 2);
	}
	return sqrt(dist);
}

//Returns a vector containing the mean of each variable in a matrix.
gsl_vector *getMeans(const gsl_matrix *mat) {
	gsl_vector *means = gsl_vector_alloc(mat->size2);
	
	for(size_t var = 0; var < mat->size2; var++) {
		double mean = 0.0;
		for(size_t point = 0; point < mat->size1; point++) {
			mean += gsl_matrix_get(mat, point, var);
		}
		mean /= mat->size1;
		gsl_vector_set(means, var, mean);
	}
	
	return means;
}

//Subtract the column means to centre the data on the origin.
gsl_matrix *centreData(const gsl_matrix *raw) {
	gsl_matrix *ret = gsl_matrix_alloc(raw->size1, raw->size2);
	
	gsl_vector *means = getMeans(raw);
	
	//Calculate the deviations by subtracting the column means from every element.
	for(size_t var = 0; var < raw->size2; var++) {
		for(size_t point = 0; point < raw->size1; point++) {
			double deviation = gsl_matrix_get(raw, point, var)
				- gsl_vector_get(means, var);
			gsl_matrix_set(ret, point, var, deviation);
		}
	}
	
	gsl_vector_free(means);
	
	return ret;
}

//Get the covariance matrix of some centred data.
gsl_matrix *getCov(const gsl_matrix *centred) {
	gsl_matrix *cov = gsl_matrix_alloc(centred->size2, centred->size2);
	
	int n = centred->size1;//The number of points.
	
	//The covariance matrix is equal to the conjugate-transpose of the centred matrix, times the centred matrix, divided by (n-1).
	//The conjugate-transpose is of course just the transpose on this real matrix.
	gsl_blas_dgemm(CblasConjTrans, CblasNoTrans, 1/((double)(n-1)), centred, centred, 0.0, cov);
	
	return cov;
}

//Get the PCA transformed and normalised data, in the fashion of Maesschalck et al. (2000).
//Requires centred data and a covariance matrix.
gsl_matrix *pcaTransform(const gsl_matrix *centred, const gsl_matrix *cov) {
	//Find the eigenthings of the covariance matrix.
	Eigen *eig = getEigen(cov);
	
	//Allocate and calculate the transformed (but unnormalised) data.
	gsl_matrix *transformed = gsl_matrix_alloc(centred->size1, centred->size2);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, centred, eig->vectors, 0.0, transformed);
	
	//Normalise the data so that each axis has unit variance.
	// Var(aX)=a^2Var(X)
	//We want Var(aX)=1
	// 1=a^2Var(X)
	// a=1/sqrt(Var(X))
	//The corresponding eigenvalue is the variance along each axis.
	//So divide by the square root of the eigenvalue.
	for(size_t var = 0; var < transformed->size2; var++) {
		for(size_t point = 0; point < transformed->size1; point++) {
			double normed;
			//Ensure that we are not dividing by zero.
			//Or trying to take the square root of a negative number.
			//Because a covariance matrix is positive semi-definite,
			//it can't have negative values.
			//So any negative values are rounding errors and should be 0.
			//TODO: Check this; I think positive semi-definite matrices may allow negative values outside the main diagonal.
			if(gsl_vector_get(eig->values, var) <= 0.0) {
				normed = 0.0;
			} else {
				normed = gsl_matrix_get(transformed, point, var)
					/ sqrt(gsl_vector_get(eig->values, var));
			}
			gsl_matrix_set(transformed, point, var, normed);
		}
	}
	
	/*for(int i = 0; i < eig->values->size; i++) {//Print the eigenvalues.
		printf("Eig: %lf\n", gsl_vector_get(eig->values, i));
	}*/
	
	freeEigen(eig);
	return transformed;
}

//Perform Hotelling's two-sample T^2-test to see if we should reject the assumption
//that the two sets of roles are drawn from the same distribution.
//Takes a vector of two matrices with the two samples, which must have been
//pca-transformed and normalised into the same co-ordinate space.
//Prints its result, and also returns the p-value.
double hotelling(vector<const gsl_matrix*> mats) {
	if(mats.size() != 2) {
		//Hotelling's is a two-sample test, but this should never trigger.
		return -1.0;
	}
	for(size_t i = 0; i < mats.size()-1; i++) {
		if(mats[i]->size2 != mats[i+1]->size2) {
			fprintf(stderr, "Hotelling's test requires that all points have the same number of dimensions.\n");
			exit(1);
		}
	}
	
	
	int p = mats[0]->size2;
	
	int n = 0;
	vector<int> ns;
	vector<gsl_vector*> means;
	for(size_t i = 0; i < mats.size(); i++) {
		means.push_back(getMeans(mats[i]));
		ns.push_back(mats[i]->size1);
		n += ns[i];
	}
	
	//TODO: Ignore any variables with mutual zeroes?
	if(n-p-1 < 0) {
		printf("Hotelling's test requires that there are more nodes in the two networks combined than components in the roles.\n");
		return 1.0;
	}
	
	//Get the Mahalanobis distance between the centres of the two distributions.
	double dist = getDist(means[0], means[1]);
	//Hotelling's t^2 statistic.
	double tsquare = dist * ((double)(ns[0] * ns[1]) / n);
	//An f statistic
	double f = tsquare * (double)(n-p-1) / ((n-2) * p);
	//TODO: Test that n-p-1 >= 0.
	
	double pvalue = gsl_cdf_fdist_Q(f, p, n-p-1);
	
	printf("Hotelling's test p-value: %lf\n", pvalue);
	
	for(size_t i = 0; i < mats.size(); i++) {
		gsl_vector_free(means[i]);
	}
	
	return pvalue;
}

//Take a vector of pointers to networks.
//Every role must have the same number of dimensions.
//Perform Principal Coordinate Analysis upon all the roles in the given networks,
//and normalise the data to give unit variance in all dimensions.
//Assume that the roles of every network are drawn from the same distribution.
//This is necessary to leave both in the same co-ordinate space after PCA.
void pca_norm_roles(vector<Network*> nets) {
	size_t num_points = 0;
	for(size_t i = 0; i < nets.size(); i++) {
		num_points += nets[i]->roles.size();
	}
	if(num_points == 0) return;
	size_t num_dims = nets[0]->roles[0].f.size();
	
	//Move the roles into a gsl matrix.
	gsl_matrix *raw = gsl_matrix_alloc(num_points, num_dims);
	size_t point = 0;
	for(size_t net = 0; net < nets.size(); net++) {
		for(size_t role = 0; role < nets[net]->roles.size(); role++) {
			for(size_t dim = 0; dim < num_dims; dim++) {
				gsl_matrix_set(raw, point, dim, nets[net]->roles[role].f[dim].frequency);
			}
			point++;
		}
	}
	
	//Centre the data and find the covariance matrix.
	gsl_matrix *centred = centreData(raw);
	gsl_matrix *cov = getCov(centred);
	
	//Perform the transformation.
	gsl_matrix *pca = pcaTransform(centred, cov);
	
	//Before we stick everything back in the roles,
	//perform Hotelling's test of the assumption.
	gsl_matrix_view *pcaMatViews = (gsl_matrix_view*)malloc(nets.size() * sizeof(gsl_matrix_view));
	vector<const gsl_matrix*> pcaMats;
	size_t row = 0;
	for(size_t i = 0; i < nets.size(); i++) {
		pcaMatViews[i] = gsl_matrix_submatrix(pca, row, 0, nets[i]->roles.size(), num_dims);
		pcaMats.push_back(&pcaMatViews[i].matrix);
		row += nets[i]->roles.size();
	}
	hotelling(pcaMats);
	fflush(stdout);//Make sure the test result gets out ASAP.
	free(pcaMatViews);
	
	//Stick everything back in the networks' roles.
	point = 0;
	for(size_t net = 0; net < nets.size(); net++) {
		for(size_t role = 0; role < nets[net]->roles.size(); role++) {
			for(size_t dim = 0; dim < num_dims; dim++) {
				nets[net]->roles[role].f[dim].frequency = gsl_matrix_get(pca, point, dim);
			}
			point++;
		}
	}
	
	//Free everything.
	gsl_matrix_free(pca);
	gsl_matrix_free(cov);
	gsl_matrix_free(centred);
	gsl_matrix_free(raw);
}
