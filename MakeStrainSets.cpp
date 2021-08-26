//Command Line Arguments:
// -i(nitial poscar location)
// -d(elta)
// -n(umber of deltas)
// -g(roup id of crystal)

//Given an initial POSCAR, creates a set of directories
//With the following naming convention:
// `pwd`/ohessSetNum/deltaNum/

//ohessSetNum depends on the type of set asked 4 (e.g. 3 for orthorhombic)
//deltaNum depends on how many deltas are asked for (e.g. nDeltas - 1)
//the - 1 is because zero strain is stored in `pwd`/init


#include <string>
#include <iomanip>
#include "Deform.h"
#include "OHESSDatabase.h"
#include "PoscarInfo.h"
#include "CrossPlatformUtils.h"

#include <iostream>

bool AlmEqu(double, double);
std::string GetPathString(int, int, double);

std::string INIT_POS_LOC = "POSCAR";
double DELTA_MAX = 0.01; //max strain factor
int N_DELTAS = 5; 
int CRYST_GROUP_ID = 8;

int main(int argc, char *argv[]){
	//CLIs
	for(int i = 0; i < argc; i++){
		if(argv[i][0] == '-' && (argv[i][1] == 'i' || argv[i][1] == 'I'))
			INIT_POS_LOC = argv[i+1];
		if (argv[i][0] == '-' && (argv[i][1] == 'd' || argv[i][1] == 'D'))
			DELTA_MAX = std::stof(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'n' || argv[i][1] == 'N'))
			N_DELTAS = std::stoi(argv[i + 1]);
		if (argv[i][0] == '-' && (argv[i][1] == 'g' || argv[i][1] == 'G'))
			CRYST_GROUP_ID = std::stoi(argv[i + 1]);
	}
	///N_DELTAS must be odd so that zero strain is included
	N_DELTAS = (N_DELTAS%2 == 0)? N_DELTAS + 1 : N_DELTAS;


	//Set up array of deltas to apply to each strain matrix
	double *deltas = (double*)malloc(N_DELTAS * sizeof(double));
	Linspace(deltas, N_DELTAS, -DELTA_MAX, +DELTA_MAX);
	
	//Read in initial POSCAR, set up a calculation for it
	Poscar initPoscar("readAll", INIT_POS_LOC.c_str());
	std::string cwd = GetCwd();
	MakePath(cwd + dl + std::to_string(CRYST_GROUP_ID) + "_-1_" + "0.00000");
	initPoscar.write(cwd + dl + std::to_string(CRYST_GROUP_ID) + "_-1_" + "0.00000" + dl + "POSCAR");

	//Loop through all of this crystal group's matrices...
	for(int i = 0; i < NumNonzeroMatrices(CRYST_GROUP_ID); i++){
		double theseEpsilons[6] = {0.0};
		//And various deltas (excluding 0.0)...
		for(int j = 0; j < N_DELTAS; j++){
			if(AlmEqu(deltas[j], 0.0))
				continue;

			///Apply delta, create new directory
			GetEpsilons(theseEpsilons, deltas[j], CRYST_GROUP_ID, i);
			std::string thisPath = cwd + dl + GetPathString(CRYST_GROUP_ID, i, deltas[j]);
			MakePath(thisPath);

			///Apply deformation, write new poscar 
			Poscar thisPoscar = initPoscar;
			DeformPoscar(thisPoscar, theseEpsilons);
			thisPoscar.write(thisPath + dl + "POSCAR");
		}
	}





	return 0;

}

//Avoids errors with floatin point arith
bool AlmEqu(double a, double b){
	return std::abs(a - b) < DELTA_MAX / 100000.0;
}

//Rules for creating new directory path
std::string GetPathString(int crystId, int matrixId, double delta){
	std::stringstream deltaStream;
	deltaStream << std::setprecision(5) << std::fixed << delta;
	std::string deltaStr = deltaStream.str();
	if(deltaStr[0] == '-')
		deltaStr[0] = 'n';
	else
		deltaStr = 'p' + deltaStr;

	return *new std::string(std::to_string(crystId) + "_" + std::to_string(matrixId) + "_" + std::string(deltaStr));
}