#pragma once
/*
*  Header file for applying deformation matrices
*  Uses instances of Poscar from PoscarInfo.h
*/

#include "PoscarInfo.h"


//Kron. Delta function: 1 if i = j else 0
double KD(int i, int j){
	return (i == j)? 1.0 : 0.0;
}

//Fills an array of size n with n equally spaced decimals from (and including) a -> b w/ b > a
void Linspace(double *arr, int n, double a, double b){
	double spacing = (b - a)/double(n - 1);
	for (int i = 0; i < n; i++){
		arr[i] = a + spacing*double(i);
	}

	return;
}

//Fills a 3x3 deformation matrix 
//ep1 - ep6 (stored in eps) are the epsilons of the strain matrix following voigt notation for the indices
//The constructed matrix follows the engineering formalism (off-diagonals are divided by 2)
//It takes the form def[3][3] = {{e1,   e6/2, e5/2}, 
//								 {e6/2, e2,   e4/2},  +   3x3 identity matrix 
//								 {e5/2, e4/2, e3}}
//therefore def[i][j] = the i'th row and the j'th column 
int EPS_SIZE = 6;
void GetDeformMatrix(double **def, const double *eps, int size=EPS_SIZE){
	///diagonals
	def[0][0] = eps[0] + 1.0;
	def[1][1] = eps[1] + 1.0;
	def[2][2] = eps[2] + 1.0;

	///upper triang
	def[0][1] = eps[5]/2.0;
	def[0][2] = eps[4]/2.0;
	def[1][2] = eps[3]/2.0;

	///symmetric lower triang
	def[1][0] = def[0][1];
	def[2][0] = def[0][2];
	def[2][1] = def[1][2];

	return;
}

//Edits a length 3 lattice vector by applying a deformation to it
//If D is the deformation matrix, then the new lattice vector a' is found by D dot a (a is the original column lattice vector)
//or a dot D (where a is the original row lattice vector).  
//Follows the notation of 'Strain and Stress: Derivation, Implementation, and Application to Organic Crystals' - Franz Knuth, 2015
//and 'High-efficiency calculation of elastic constants enhanced by the optimized strain - matrix sets' - Zhong-Li Liu, 2020
void DeformVec(double *vec, double **defMat){
	const double vecOrig[3] = {vec[0], vec[1], vec[2]};
	
	for(int i = 0; i < 3; i++){ ///i = {1, 2, 3} -> {x, y, z} (NOT voigt notation here)
		vec[i] = 0.0;
		for(int j = 0; j < 3; j++){ ///ditto above comment
			vec[i] += defMat[i][j] * vecOrig[j];
		}
	}
	
	return;
}

//Applies deformation matrix to this poscar instance
//Needs a length 6 array of epsilons following voigt convention (minding the 0-indexing of C++ i.e. epsilons[3] -> voigt epsilon 4)
void DeformPoscar(Poscar &pos, const double *epsilons){
	pos.applyUsc();
	pos.convertToDirect();

	//Initialize deformation matrix (with malloc for some reason?  I can't remember why I did this but I'm not messing with it now)
	double **deformMatrix = (double**)malloc(3*sizeof(double*));
	for(int i = 0; i < 3; i++)
		deformMatrix[i] = (double*)malloc(3*sizeof(double));
	GetDeformMatrix(deformMatrix, epsilons);

	//Edit lattice vectors
	double newVecA[3] = {pos.superCellVectorA[0], pos.superCellVectorA[1], pos.superCellVectorA[2]};
	DeformVec(newVecA, deformMatrix);
	double newVecB[3] = {pos.superCellVectorB[0], pos.superCellVectorB[1], pos.superCellVectorB[2]};
	DeformVec(newVecB, deformMatrix);
	double newVecC[3] = {pos.superCellVectorC[0], pos.superCellVectorC[1], pos.superCellVectorC[2]};
	DeformVec(newVecC, deformMatrix);
	for(int i = 0; i < 3; i++){
		pos.superCellVectorA[i] = newVecA[i];
		pos.superCellVectorB[i] = newVecB[i];
		pos.superCellVectorC[i] = newVecC[i];
	}

	return;
}



