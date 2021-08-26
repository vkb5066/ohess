#pragma once
/*
* Header file contiaining all strain sets necessary for ohess
* OHESS defined in 'High-efficiency calculation of elastic constants enhanced by the optimized strain - matrix sets' - Zhong-Li Liu, 2020
*/ 

//Format is the following:
// OHESS_DB[n] = LAT_SYS (0 <= n < 9)
// LAT_SYS[m] = EPS (0 <= m < 6)
// EPS[o] = 1.0 if nonzero or zero if zero (0 <= o < 6)
//   i.e. multiplynig a delta to the fators in a matrix is perfectly fine

//Database Keys:
/*
* 0 = Cubic
* 1 = Hexagonal
* 2 = Rhomb I
* 3 = Rhomb II
* 4 = Tetrag I
* 5 = Tetrag II
* 6 = Ortho
* 7 = Mono
* 8 = Tri
*/ 

extern const double OHESS_DB[9][6][6] = 
{	//---CUBIC--- (1 nonzero)
	{	
		{1.0, 0.0, 0.0, 1.0, 0.0, 0.0}, //c11, c12, c44
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	},
	//---HEXAG--- (2 nonzeros)
	{	
		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c11, c12, c13
		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c44
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	},
	//---RHOMB I--- (2 nonzeros)
	{	
		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c1i, i = {1-4}
		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c44
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	},
	//---RHOMB II--- (2 nonzeros)
	{	
		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c1i, i = {1-5}
		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c44
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	},
	//---TETRAG I--- (2 nonzeros)
	{	
		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c11, c12, c13
		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c44
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	},
	//---TETRAG II--- (2 nonzeros)
	{	
		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c1i, i = {1-3, 6}
		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c44
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	},
	//---ORTHO--- (3 nonzeros)
	{	
		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c11, c12, c13
		{0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, //c22, c23
		{0.0, 0.0, 1.0, 1.0, 1.0, 1.0}, //cii, i = {3-6}
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	},
	//---MONO--- (4 nonzeros)
	{	
		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c1i, i = {1-3, 6}
		{0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, //c22, c23, c26
		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c36, c44, c45
		{0.0, 0.0, 0.0, 0.0, 1.0, 1.0}, //c55, c66
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	},
	//---TRICLIN--- (6 nonzeros)
	{	
		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c1i, i = {1-6}
		{0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, //c2i, i = {2-6} 
		{0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, //c3i, i = {3-6} 
		{0.0, 0.0, 0.0, 1.0, 0.0, 0.0}, //c44, c45, c46 
		{0.0, 0.0, 0.0, 0.0, 1.0, 0.0}, //c55, c56 
		{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}  //c66
	}
};

//Returns the number of nonzero epsilon matrices for a given lattice system
//i.e. the argument 6 returns 3 since the orthorhombic system has 3 sets of epsilons
int NumNonzeroMatrices(int latSysId){
	switch(latSysId){
		case 0:
			return 1;
		case 1:
			return 2;
		case 2:
			return 2;
		case 3:
			return 2;
		case 4:
			return 2;
		case 5:
			return 2;
		case 6:
			return 3;
		case 7:
			return 4;
		case 8:
			return 6;
		default:
			return -1;
	}
}

//Copies a given strain matrix to an array
//Applies a multiplication factor (should be on the order of +/-0.01 to +/-0.001) to prepare to apply strains
void GetEpsilons(double *arrToModify, double multFac, int latId, int epsId){
	for(int i = 0; i < 6; i++){
		arrToModify[i] = multFac * OHESS_DB[latId][epsId][i];
	}
}





