#pragma once
//Header file for reading in POSCAR and CONTCAR files.  
//As of right now, file formats that include atomic velocities at the bottom are not supported, and will likley cause errors
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <sstream>

//Globals:  stuff to change paramaters as needed
extern double GENERAL_BOND_DISTANCE = 3.5; ///(angstroms) distance between any two atoms for them to be considered bonded (works well for POSCARS with only a few atom types)
//extern std::string PATH_TO_FILES = "/storage/work/vkb5066/scriptTests/tst19/";
extern std::string PATH_TO_FILES = "C://Users//baron//Desktop//";
extern std::string DEFAULT_POSCAR_OPEN_PATH = (PATH_TO_FILES + "POSCAR").c_str();
extern std::string DEFAULT_BONDINFO_OPEN_PATH = (PATH_TO_FILES + "bondDataInfo").c_str();


//Main variable type for this class.  Holds information about individual atoms
struct Coords
{
	std::string atomType; ///holds the type of atom being described
	long double a, b, c; ///atomic coordinates 
	std::string flags; ///string for the optional molecular dynamics flags after the coordinates.  Can also hold a blank string to avoid annoying errors
	int id; ///useful way for IDing atoms in the case that an iterative loop is not enough

	std::string extraInfo; ///string for tags that may be useful for other calculations.  Not used in this header file
	bool extraVal; ///bool that may be useful for other calculations

	//Default constructor
	Coords()
	{
		atomType = "undf"; ///can probably be left alone, but may be useful
		a = -99;
		b = -99;
		c = -99;
		flags = " Uninitialized instance...you should not see this. : Coords - Default Constructor";
		id = -1;
		extraInfo = "undf : Coords: - Default Constructor";
		extraVal = false;
	}

	//Param constructor
	Coords (long double a_, long double b_, long double c_)
	{
		Coords();
		a = a_;
		b = b_;
		c = c_;
	}

	//Operator declerations
	Coords& operator= (const Coords &coords);
	bool operator== (const Coords &) const;
	bool operator&= (const Coords &) const;
	bool operator!= (const Coords &) const;
};
double dist_(const Coords, const Coords);
Coords& Coords::operator= (const Coords &coords)
{
	atomType = coords.atomType;
	a = coords.a;
	b = coords.b;
	c = coords.c;
	flags = coords.flags;
	id = coords.id;
	extraInfo = coords.extraInfo;
	return *this;
}
bool Coords::operator ==(const Coords &coords) const ///Compare by distance between atoms
{
	return (dist_(*this, coords) < 0.01); ///TODO:  try to find a better way to determine if two atoms are close to eachother
}
bool Coords::operator &=(const Coords &coords) const ///Compare by ID
{
	return id == coords.id;
}
bool Coords::operator !=(const Coords & coords) const
{
	return !(*this == coords);
}

//Extra variable type.  Holds information about pairs (not necessairily bonds) of atoms.  Max supported is two atoms per pair
struct atomPair
{
	Coords pairedAtoms[2]; ///info about the two related atoms
	double lenBetween; ///length between the two paired atoms
	double genNum; ///general number for various calculations.  should be reset at the end of every function, since it is used for multiple things
	std::string type; ///type of pair, not type of atom. for example, two atoms A and B should create a type of AB
	std::string type2; ///type of pair, but the atoms are flipped around
	std::string extraInfo;

	//Default class constructor.  Not overly useful, but may be needed in some implementations
	atomPair()
	{
		Coords a, b; ///leave these uninitialized -- Poscar class constructor, or a call to the initialization, should set values automatically
		pairedAtoms[0] = a;
		pairedAtoms[1] = b;

		genNum = 0;
		lenBetween = 0;
		type = "uninitialized...you should not see this. : Atompair - Default Constructor";
		type2 = "uninitialized...you should not see this. : Atompair - Default Constructor";
		extraInfo = "undf : AtomPair - Default Constructor";
	}

	//Paramaterized class constructor.  Uses the info from the two input atoms, ***under the assumption that they are initialized already***
	atomPair(const Coords atomA, const Coords atomB)
	{
		pairedAtoms[0] = atomA;
		pairedAtoms[1] = atomB;

		genNum = 0;
		lenBetween = dist_(atomA, atomB);
		type = (atomA.atomType + atomB.atomType).c_str(); ///should add a sort by alphabetical order so C-H bonds and H-C bonds are considered the same thing
		type2 = (atomB.atomType + atomA.atomType).c_str();
		extraInfo = "undf";
	}

	//Operator declerations
	atomPair& operator= (const atomPair &pair);
	bool operator== (const atomPair &) const;
	bool operator&= (const atomPair &) const;
	bool operator^= (const atomPair &) const;
	bool operator!= (const atomPair &) const;
};
atomPair& atomPair::operator= (const atomPair &pair)
{
	pairedAtoms[0] = pair.pairedAtoms[0];
	pairedAtoms[1] = pair.pairedAtoms[1];
	lenBetween = pair.lenBetween;
	genNum = pair.genNum;
	type = pair.type;
	type2 = pair.type2;
	extraInfo = pair.extraInfo;
	return *this;
}
bool atomPair::operator ==(const atomPair &atompair) const ///Compare by Position
{
	return (((pairedAtoms[0] == atompair.pairedAtoms[0]) && (pairedAtoms[1] == atompair.pairedAtoms[1])) || ((pairedAtoms[0] == atompair.pairedAtoms[1]) && (pairedAtoms[1] == atompair.pairedAtoms[0])));
}
bool atomPair::operator &=(const atomPair &atompair) const ///Compare by ID
{
	return (((pairedAtoms[0].id == atompair.pairedAtoms[0].id) && (pairedAtoms[1].id == atompair.pairedAtoms[1].id)) || ((pairedAtoms[0].id == atompair.pairedAtoms[1].id) && (pairedAtoms[1].id == atompair.pairedAtoms[0].id)));
}
bool atomPair::operator ^=(const atomPair &atompair) const ///Compare by type
{
	return((type == atompair.type) || (type2 == atompair.type2) || (type == atompair.type2) || (type2 == atompair.type));
}
bool atomPair::operator !=(const atomPair & atompair) const
{
	return !(*this == atompair);
}

struct bondInfo
{
	std::string type1;
	std::string type2;
	double physBondDist;

	bondInfo()
	{
		type1 = "bondInfo(): undf";
		type2 = "bondInfo(): undf";
		physBondDist = -99;
	}

	bondInfo(std::string one, std::string two)
	{
		type1 = (one + two).c_str();
		type2 = (two + one).c_str();
		physBondDist = -99;
	}

	bondInfo(std::string one, std::string two, double dist)
	{
		type1 = (one + two).c_str();
		type2 = (two + one).c_str();
		physBondDist = dist;
	}
};

//Comparison rule for sorting strings
bool stringSort_(std::string a, std::string b)
{
	return a < b;
}

//-------------------------------------------------------------------------------------------------------------------

//Class Decleration
class Poscar
{
	public:
		//Variables used in this header file
		std::string defaultPath = DEFAULT_POSCAR_OPEN_PATH.c_str();
		std::string defaultWritePath = (DEFAULT_POSCAR_OPEN_PATH + "write").c_str();
		std::string infilePath;
		std::string fileTitle; ///first line of file that is ignored by VASP
		std::string modelType; ///defines if model is bulk or molecular
		double universalScaleFactor; ///all atomic coordinates and lattice vectors are scaled by this.  A negative value is read as a predefined cell volume
		long double superCellVectorA [3];  ///-----      
		long double superCellVectorB [3];  ///      > defines the supercell
		long double superCellVectorC [3];  ///-----
		double volume; ///in angstroms^3
		std::vector <std::string> atomTypes; ///vector to hold arbitrary number of atom types ('C', 'Ag', etc.)
		std::vector <int> atomTypeNums; ///vector (parllel to above) to hold the number of atoms of each type
		bool selectiveDynamicsTag; ///false if the selective dynamics tag is not there
		bool directTag = false, cartesianTag = false; ///true if the system is in direct or cartesian coordinates, respectivly
		std::vector <Coords> atomCoords; ///vector of all Coords (defined before the class)
		std::vector <atomPair> atomPairs; ///vector of every atom's relation to every other atom.  get rid of this if the code takes too long to run
		std::vector <atomPair> atomBonds; ///similar to atomPairs, but with a distance restriction on what atoms are related to one another
		std::vector <bondInfo> bondInfoVect;

		//Variables used in 'Calculations' header file
		std::string bondDataInfoLoc = DEFAULT_BONDINFO_OPEN_PATH.c_str();

		//Constructor Declerations
		Poscar();
		Poscar(std::string);
		Poscar(std::string, std::string);

		//Functions Declerations included in this header file
		void fetchFileTitle();
		void fetchModelType();
		void fetchUniversalScaleFactor();
		void fetchSuperCellVectors();
		void fetchSupercellVolume();
		void fetchAtomTypes();
		void fetchAtomTypeNums();
		void fetchSelectiveDynamicsTag();
		void fetchDirectTag();
		void fetchCartesianTag();
		void fetchAtomCoords();
		void applyUsc(); ///apply universal scaling constant
		void convertToDirect();
		void convertToDirect(Coords &crds);
		void convertToDirect(std::vector <Coords> &crds);
		void convertToCartesian();
		void convertToCartesian(Coords &crds);
		void convertToCartesian(std::vector <Coords> &crds);
		void removeTaggedAtoms(std::string);
		void removeDuplicates(); ///note:  NOT automatically called upon construction of class instance.  need to manually call this if you want it to be used
		void removeDuplicates(double);
		void removeDuplicates(std::vector <Coords> &);
		void removeDuplicates(std::string); ///paramatarized version.  can removes by coords ("coords") or by atom ID ("id") or non-original atoms [for use with supercell extend] ("unoriginal")
		void updateAtomTypes();				///or by atoms on the outside edges / or past the supercell ("edge")
		void updateAtomCounts();
		void updateAtomCounts(std::vector <Coords>);
		void updateAtomCoords();
		void assertAtomOrder(std::vector <std::string>);
		void updateAll();
		void copyFormatting(const Poscar); ///copies the POSCAR formatting of the input file (order of atom types, order of atoms)
		void print();
		void print(std::string, std::string);
		void write();
		void write(std::string);
		void clearEmptyValues();

		Poscar& operator= (const Poscar &poscar);

		//Function Declerations included in the 'Calculations' header file
		void extendSupercell(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int, unsigned int);
		std::vector <Coords> allPeriodicImages(std::string, int);
		std::vector <Coords> allPeriodicImages(double);
		double grapheneAreaApprox(); ///assuming the POSCAR describes a 4-vertex plane, calculates the surface area of the carbons
		void translateAtoms(double, double, double);
		void applyStrain(bool, bool, bool, int, double); ///strains models
		void fetchAtomPairs(double, int); ///NOT automatically called in class constructor
		void fetchAtomBonds(); ///NOT automatically called in class constructor
		void fetchAtomBonds(std::string);
		void fetchBondInfoVect(std::string);
		void restorePositions(); ///Puts any atoms outside of the supercell back into the supercell
		//void writeBondInfo(std::string);
};

//Class Constructors-----------------------------------------------------------------------
//Default class constructor; user will have to supply most of the info somehow
Poscar::Poscar()
{
	infilePath = defaultPath;
	std::string str = "nothing read in: Poscar - Default Constructor";
	fileTitle = str;
	modelType = "undf";
	universalScaleFactor = 0;
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = -1;
		superCellVectorB[i] = -1;
		superCellVectorC[i] = -1;
	}
	volume = -1;
	atomTypes.push_back(str);
	atomTypeNums.push_back(-1);
	selectiveDynamicsTag = 0;
	directTag = 0;
	cartesianTag = 0;
	atomCoords;
	atomPairs;
	atomBonds;
}

//Default class constructor; user will have to supply most of the info somehow; but with the filepath defined
Poscar::Poscar(std::string path)
{
	infilePath = path;
	defaultWritePath = (infilePath + "write").c_str();

	std::string str = "nothing read in: Poscar - string Constructor";
	fileTitle = str;
	modelType = "undf";
	universalScaleFactor = 0;
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = -1;
		superCellVectorB[i] = -1;
		superCellVectorC[i] = -1;
	}
	volume = -1;
	atomTypes.push_back(str);
	atomTypeNums.push_back(-1);
	selectiveDynamicsTag = 0;
	directTag = 0;
	cartesianTag = 0;
	atomCoords;
	atomPairs;
	atomBonds;
}

//Special (useful) class constructor.  Input: (string) instructions["readHead"] or ["readAll"], (string) file path
Poscar::Poscar(std::string instructions, std::string path)
{
	//Reads only the head of the file.  Does NOT return an iterator to the beginning of atomic coordinates; user will have to do this before reading those in
	if (instructions == "ReadHead" || instructions == "readHead")
	{
		infilePath = path;
		defaultWritePath = (infilePath + "write").c_str();

		fetchFileTitle();
		fetchUniversalScaleFactor();
		fetchSuperCellVectors();
		fetchSupercellVolume();
		fetchAtomTypes();
		fetchModelType();
		fetchAtomTypeNums();
		fetchSelectiveDynamicsTag();
		fetchDirectTag();
		fetchCartesianTag();
	}
	else
	//Reads the entire file, storing atomic coordinates in a vector
	if (instructions == "ReadAll" || instructions == "readAll")
	{
		infilePath = path;
		defaultWritePath = (infilePath + "write").c_str();

		fetchFileTitle();
		fetchUniversalScaleFactor();
		fetchSuperCellVectors();
		fetchSupercellVolume();
		fetchAtomTypes();
		fetchModelType();
		fetchAtomTypeNums();
		fetchSelectiveDynamicsTag();
		fetchDirectTag();
		fetchCartesianTag();
		fetchAtomCoords();
		atomPairs;
		atomBonds;
	}

	else
		std::cout << "Poscar Constructor:  unknown instructions argument\n";
}


//Destructor (not needed unless dynamic memory allocation is used.  As of right now, it isnt.)
//~Poscar() {}

//-----------------------------------------------------------------------------------------
//Overloaded = sign.  This "should" be very useful
Poscar& Poscar::operator= (const Poscar &poscar)
{
	infilePath = poscar.infilePath;
	defaultWritePath = poscar.defaultWritePath;
	fileTitle = poscar.fileTitle;
	universalScaleFactor = poscar.universalScaleFactor;
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = poscar.superCellVectorA[i];
		superCellVectorB[i] = poscar.superCellVectorB[i];
		superCellVectorC[i] = poscar.superCellVectorC[i];
	}
	volume = poscar.volume;
	atomTypes = poscar.atomTypes;
	atomTypeNums = poscar.atomTypeNums;
	modelType = poscar.modelType;
	selectiveDynamicsTag = poscar.selectiveDynamicsTag;
	directTag = poscar.directTag;
	cartesianTag = poscar.cartesianTag;
	atomCoords = poscar.atomCoords;
	atomPairs = poscar.atomPairs;
	atomBonds = poscar.atomBonds;
	return *this;
}

//----------------------------------------------------------------------------------------------------------------------

//Non-Member Function Definitions----------------------------------------------------------------------------------------

std::string readNthLine(const std::string filename, int n)
{
		std::ifstream infile(filename.c_str());
		std::string s;

		//Make this faster
		s.reserve(1024);

		//skip n lines
		for (int i = 0; i < n; ++i)
			std::getline(infile, s);

		std::getline(infile, s);
		return s;
}

double angleBetween (std::vector <double> a, std::vector <double> b)
{
	return acos(((a[0]*b[0]) + (a[1]*b[1]) + (a[2] * b[2]))/(sqrt((a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]))*sqrt((b[0] * b[0]) + (b[1] * b[1]) + (b[2] * b[2]))));
}

double angleBetween(long double a [3], long double b [3])
{
	return acos(((a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2])) / (sqrt((a[0] * a[0]) + (a[1] * a[1]) + (a[2] * a[2]))*sqrt((b[0] * b[0]) + (b[1] * b[1]) + (b[2] * b[2]))));
}

//Distance between two points in R3.  Why can I not use the one defined in Calculations.h?  Even though i'm trying to include it?  I dont know
double dist_(const Coords coordA, const Coords coordB)
{
	return sqrt((((coordB.a - coordA.a)*(coordB.a - coordA.a))) + (((coordB.b - coordA.b)*(coordB.b - coordA.b))) + (((coordB.c - coordA.c)*(coordB.c - coordA.c))));
}

//Formula taken from wikipedia: https://en.wikipedia.org/wiki/Fractional_coordinates
void setTransformMatrixToDirect (std::vector <double> vectA, std::vector <double> vectB, std::vector <double> vectC, std::vector <double> &toReturn)
{
	toReturn.clear();
	//syntax for vector positions in reference to the matrix:
	//        | [0]  [1]  [2] |
	//        | [3]  [4]  [5] |
	//        | [6]  [7]  [8] |

	//Definitions (if you are debugging this: 1) i'm sorry.  2) you should look at the wiki link for the matrix formula used there)
	double a = sqrt((vectA[0] * vectA[0]) + (vectA[1] * vectA[1]) + (vectA[2] * vectA[2]));
	double b = sqrt((vectB[0] * vectB[0]) + (vectB[1] * vectB[1]) + (vectB[2] * vectB[2]));
	double c = sqrt((vectC[0] * vectC[0]) + (vectC[1] * vectC[1]) + (vectC[2] * vectC[2]));
	double alpha = angleBetween(vectB, vectC); //the angle between b and c
	double beta = angleBetween(vectA, vectC); //the angle between a and c
	double gamma = angleBetween(vectA, vectB); //the angle between a and b
	double volume = a*b*c*sqrt(1 - (cos(alpha) * cos(alpha)) - (cos(beta) * cos(beta)) - (cos(gamma) * cos(gamma)) + 2*cos(alpha)*cos(beta)*cos(gamma));
	
	//Filling vector / matrix
	toReturn.push_back(1 / a); //[0]
	toReturn.push_back(-cos(gamma) / (a * sin(gamma))); //[1]
	toReturn.push_back(b * c * ((cos(alpha) * cos(gamma) - cos(beta))/(volume * sin(gamma)))); //[2]
	toReturn.push_back(0); //[3]
	toReturn.push_back(1 / (b * sin(gamma))); //[4]
	toReturn.push_back(a * c * ((cos(beta) * cos(gamma) - cos(alpha)) / (volume * sin(gamma)))); //[5]
	toReturn.push_back(0); //[6]
	toReturn.push_back(0); //[7]
	toReturn.push_back(a * b * sin(gamma)/(volume)); //[8]
}

//Taken from https://stackoverflow.com/questions/8425214/splitting-string-into-a-vectorstring-of-words
std::vector<std::string> split(const std::string& s)
{
	std::vector<std::string> ret;
	typedef std::string::size_type string_size;
	string_size i = 0;

	// invariant: we have processed characters [original value of i, i) 
	while (i != s.size()) {
		// ignore leading blanks
		// invariant: characters in range [original i, current i) are all spaces
		while (i != s.size() && isspace(s[i]))
			++i;

		// find end of next word
		string_size j = i;
		// invariant: none of the characters in range [original j, current j)is a space
		while (j != s.size() && !isspace(s[j]))
			j++;

		// if we found some nonwhitespace characters 
		if (i != j) {
			// copy from s starting at i and taking j - i chars
			ret.push_back(s.substr(i, j - i));
			i = j;
		}
	}
	return ret;
}

bool empty(const std::string s)
{
	std::string s_  = s.c_str(); ///c_str() because unix dosent like "normal" strings
	if (s_ == "" || s_ == " " || s_ == "\t" || s_.length() == 0)
		return 1;
	else return 0;
}


//Member Function Definitions:-------------------------------------------------------------------------------------------------------------------

void Poscar::fetchFileTitle()
{
	//Open file:
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "fetchFileTitle(): File could not be found...expect errors\n";

	//Get first line, update the class title
	getline(infile, fileTitle);
	fileTitle = fileTitle.c_str(); ///c_str() because unix dosent like "normal" strings
	
	infile.close();
}

//Warning: works based off of whether the input file has carbons (molecular model) or no carbons (bulk model).  Don't rely on this to work in general.
void Poscar::fetchModelType()
{
	if (atomTypes.size() <= 1)
		modelType = "Bulk";
	else
		modelType = "Molecular";
}

void Poscar::fetchUniversalScaleFactor()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "fetchUnivScaleFactor(): File could not be found...expect errors\n";
	
	//move to line before line to read in
	int n = 1; ///number of lines to skip
	for (int i = 0; i < n; i++)
	{
		for(;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update class param
	infile >> universalScaleFactor;

	infile.close();
}

void Poscar::fetchSuperCellVectors()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "fetchSupercellVectors(): File could not be found...expect errors\n";

	//move to line before line to read in
	int n = 2; ///number of lines to skip
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update class params
	for (int i = 0; i < 3; i++)
		infile >> superCellVectorA[i];
	for (int i = 0; i < 3; i++)
		infile >> superCellVectorB[i];
	for (int i = 0; i < 3; i++)
		infile >> superCellVectorC[i];

	infile.close();
}

void Poscar::fetchSupercellVolume()
{
	double a = sqrt((superCellVectorA[0] * superCellVectorA[0]) + (superCellVectorA[1] * superCellVectorA[1]) + (superCellVectorA[2] * superCellVectorA[2]));
	double b = sqrt((superCellVectorB[0] * superCellVectorB[0]) + (superCellVectorB[1] * superCellVectorB[1]) + (superCellVectorB[2] * superCellVectorB[2]));
	double c = sqrt((superCellVectorC[0] * superCellVectorC[0]) + (superCellVectorC[1] * superCellVectorC[1]) + (superCellVectorC[2] * superCellVectorC[2]));
	double alpha = angleBetween(superCellVectorB, superCellVectorC); //the angle between b and c
	double beta = angleBetween(superCellVectorA, superCellVectorC); //the angle between a and c
	double gamma = angleBetween(superCellVectorA, superCellVectorB); //the angle between a and b
	volume = a*b*c*sqrt(1 - (cos(alpha) * cos(alpha)) - (cos(beta) * cos(beta)) - (cos(gamma) * cos(gamma)) + 2 * cos(alpha)*cos(beta)*cos(gamma));
}

void Poscar::fetchAtomTypes()
{
	if (!atomTypes.empty())
		return;

	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "fetchAtomTypes(): File could not be found...expect errors\n";

	//move to line before line to read in
	int n = 5; ///number of lines to skip
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update params:
		//Get the non-blank line that holds all of the atom types
	std::string atomTypeLine;
	for (;;)
	{
		getline(infile, atomTypeLine);
		if (!empty(atomTypeLine.c_str()))
			break;
	}
	atomTypes = split(atomTypeLine);

	infile.close();
}

void Poscar::fetchAtomTypeNums()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "fetchAtomTypeNums(): File could not be found...expect errors\n";

	//move to line before line to read in
	int n = 6; ///number of lines to skip
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update params
	for (int i = 0; i < atomTypes.size(); i++)
	{	
		int n;
		infile >> n; ///should be in line with the atom types
		atomTypeNums.push_back(n);
	}

	infile.close();
}

void Poscar::fetchSelectiveDynamicsTag()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "fetchSelectiveDynamicsTag(): File could not be found...expect errors\n";

	//move to line before line to read in
	int n = 7; ///number of lines to skip
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update param
	std::string s;
	infile >> s;
	if ((s[0] == 'S' || s[0] == 's') && (s[1] == 'E' || s[1] == 'e'))
		selectiveDynamicsTag = 1;
	else selectiveDynamicsTag = 0;

	infile.close();
}

void Poscar::fetchDirectTag()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "fetchDirectTag(): File could not be found...expect errors\n";

	//need to skip different number of lines depending on whether the 'Sel' tag is there or not, since it takes up an entire line
	int n;
	if (selectiveDynamicsTag)
		 n = 8;
	else
		 n = 7;
	
	//move to line before line to read in
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update param
	std::string s;
	infile >> s;
	if ((s[0] == 'D' || s[0] == 'd') && (s[1] == 'i' || s[1] == 'i'))
	{
		directTag = 1;
		cartesianTag = 0;
	}
	else 
		directTag = 0;
	infile.close();
}

void Poscar::fetchCartesianTag()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "fetchCartesianTag(): File could not be found...expect errors\n";

	//need to skip different number of lines depending on whether the 'Sel' tag is there or not, since it takes up an entire line
	int n;
	if (selectiveDynamicsTag)
		n = 8;
	else
		n = 7;

	//move to line before line to read in
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Update param
	std::string s;
	infile >> s;
	if ((s[0] == 'C' || s[0] == 'c') && (s[1] == 'A' || s[1] == 'a'))
	{
		cartesianTag = 1;
		directTag = 0;
	}
	else cartesianTag = 0;

	infile.close();
}

void Poscar::fetchAtomCoords()
{
	//Open file
	std::ifstream infile(infilePath.c_str());
	if (infile.fail())
		std::cout << "fetchAtomCoords(): File could not be found...expect errors\n";

	//need to skip different number of lines depending on whether the 'Sel' tag is there or not, since it takes up an entire line
	int n;
	if (selectiveDynamicsTag)
		n = 9;
	else
		n = 8;

	//move to line before line to read in
	for (int i = 0; i < n; i++)
	{
		for (;;)
		{
			std::string s;
			getline(infile, s);
			if (!empty(s.c_str()))
				break;
		}
	}

	//Over-complicated reading in of atomic coordinates so that each atom can be assigned its correct atom type
	for (int i = 0; i < atomTypes.size(); i++)
	{
		for (int j = atomTypeNums[i]; j > 0; j--) ///if there will be any error in this code, it will be here, with this for loop's termination rule
		{
			Coords crds;
			infile >> crds.a >> crds.b >> crds.c;
			getline(infile, crds.flags);
			crds.atomType = atomTypes[i];
			atomCoords.push_back(crds);
		}
	}

	//And setting ID numbers for the atoms
	for (int i = 0; i < atomCoords.size(); i++)
		atomCoords[i].id = i;

	infile.close();
}

/*
void Poscar::convertToDirect()
{
	if (directTag)  ///leaving this here just in case you'd want to do something special if file is already in cartesian
	{
		return;
	} 
	else
		if (cartesianTag)
		{
			std::vector <double> transformMatrix;
			//Somehow managed to set the supercell vectors as arrays, which is annoying to work with.  Put them in vectors:
			//May be worth fixing this later to improve performance...for now, i just have the vectors in a limited scope

			//Fill transform matrix; syntax is in the above function comment
			{	
				std::vector <double> vA, vB, vC;
				for (int i = 0; i < 3; i++)
				{
					vA.push_back(universalScaleFactor * superCellVectorA[i]);
					vB.push_back(universalScaleFactor * superCellVectorB[i]);
					vC.push_back(universalScaleFactor * superCellVectorC[i]);
				}
				setTransformMatrixToDirect(vA, vB, vC, transformMatrix);
			}
			//Conversion math
			for (int i = 0; i < atomCoords.size(); i++)
			{
				atomCoords[i].a = universalScaleFactor * ((transformMatrix[0] * atomCoords[i].a) + (transformMatrix[1] * atomCoords[i].b) + (transformMatrix[2] * atomCoords[i].c));
				atomCoords[i].b = universalScaleFactor * ((transformMatrix[3] * atomCoords[i].a) + (transformMatrix[4] * atomCoords[i].b) + (transformMatrix[5] * atomCoords[i].c));
				atomCoords[i].c = universalScaleFactor * ((transformMatrix[6] * atomCoords[i].a) + (transformMatrix[7] * atomCoords[i].b) + (transformMatrix[8] * atomCoords[i].c));
			}

			//Update tags
			for (int i = 0; i < 3; i++)
			{
				superCellVectorA[i] = universalScaleFactor * superCellVectorA[i];
				superCellVectorB[i] = universalScaleFactor * superCellVectorB[i];
				superCellVectorC[i] = universalScaleFactor * superCellVectorC[i];
			}
			universalScaleFactor = 1.0;
			directTag = true;
			cartesianTag = false;
		}
}

void Poscar::convertToDirect(std::vector <Coords> &crds)
{
			std::vector <double> transformMatrix;
			//Somehow managed to set the supercell vectors as arrays, which is annoying to work with.  Put them in vectors:
			//May be worth fixing this later to improve performance...for now, i just have the vectors in a limited scope

			//Fill transform matrix; syntax is in the above function comment
			{
				std::vector <double> vA, vB, vC;
				for (int i = 0; i < 3; i++)
				{
					vA.push_back(universalScaleFactor * superCellVectorA[i]);
					vB.push_back(universalScaleFactor * superCellVectorB[i]);
					vC.push_back(universalScaleFactor * superCellVectorC[i]);
				}
				setTransformMatrixToDirect(vA, vB, vC, transformMatrix);
			}
			//Conversion math
			for (int i = 0; i < crds.size(); i++)
			{
				crds[i].a = universalScaleFactor * ((transformMatrix[0] * crds[i].a) + (transformMatrix[1] * crds[i].b) + (transformMatrix[2] * crds[i].c));
				crds[i].b = universalScaleFactor * ((transformMatrix[3] * crds[i].a) + (transformMatrix[4] * crds[i].b) + (transformMatrix[5] * crds[i].c));
				crds[i].c = universalScaleFactor * ((transformMatrix[6] * crds[i].a) + (transformMatrix[7] * crds[i].b) + (transformMatrix[8] * crds[i].c));
			}

			//Update tags
			for (int i = 0; i < 3; i++)
			{
				superCellVectorA[i] = universalScaleFactor * superCellVectorA[i];
				superCellVectorB[i] = universalScaleFactor * superCellVectorB[i];
				superCellVectorC[i] = universalScaleFactor * superCellVectorC[i];
			}
			universalScaleFactor = 1.0;
}
*/
void Poscar::applyUsc()
{
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = superCellVectorA[i] * universalScaleFactor;
		superCellVectorB[i] = superCellVectorB[i] * universalScaleFactor;
		superCellVectorC[i] = superCellVectorC[i] * universalScaleFactor;
	}
	if(cartesianTag)
	{
		for (int i = 0; i < atomCoords.size(); i++)
		{
			atomCoords[i].a = atomCoords[i].a * universalScaleFactor;
			atomCoords[i].b = atomCoords[i].b * universalScaleFactor;
			atomCoords[i].c = atomCoords[i].c * universalScaleFactor;
		}
	}
	universalScaleFactor = 1.0;
}
void Poscar::convertToDirect()
{
	if (directTag)
	{
		return;
	}
	else
		if (cartesianTag)
		{
			///apply universal scale factor to original supercell vectors, find the determinant, and transpose of the matrix made by them
			for (int i = 0; i < 3; i++)
			{
				superCellVectorA[i] = superCellVectorA[i] * universalScaleFactor;
				superCellVectorB[i] = superCellVectorB[i] * universalScaleFactor;
				superCellVectorC[i] = superCellVectorC[i] * universalScaleFactor;
			}
			for (int i = 0; i < atomCoords.size(); i++)
			{
				atomCoords[i].a = atomCoords[i].a * universalScaleFactor;
				atomCoords[i].b = atomCoords[i].b * universalScaleFactor;
				atomCoords[i].c = atomCoords[i].c * universalScaleFactor;
			}
			universalScaleFactor = 1.0;

			//Determinant:
			long double detOrig = (superCellVectorA[0] * ((superCellVectorB[1] * superCellVectorC[2]) - (superCellVectorB[2] * superCellVectorC[1])))
				- (superCellVectorB[0] * ((superCellVectorA[1] * superCellVectorC[2]) - (superCellVectorA[2] * superCellVectorC[1])))
				+ (superCellVectorC[0] * ((superCellVectorB[2] * superCellVectorA[1]) - (superCellVectorB[1] * superCellVectorA[2])));
			if (detOrig == 0)
			{
				std::cout << "converToCartesian(): error - determinant was found to be zero\n";
				return;
			}

			//Transpose of the matrix used to convert from direct to cartesian turns out to just be the matrix of vectors as written in the head of the input file:
			long double trans [9];
			for (int i = 0; i < 3; i++)
				trans[i] = superCellVectorA[i];
			for (int i = 3; i < 6; i++)
				trans[i] = superCellVectorB[i-3];
			for (int i = 6; i < 9; i++)
				trans[i] = superCellVectorC[i-6];

			//Fill the cofactor matrix with the determinants of the transposed matrix, flipping the signs when necessary
			long double coFact [9];
			coFact[0] = ((trans[4] * trans[8]) - (trans[5] * trans[7]));
			coFact[1] = -((trans[3] * trans[8]) - (trans[5] * trans[6]));
			coFact[2] = ((trans[3] * trans[7]) - (trans[4] * trans[6]));
			coFact[3] = -((trans[1] * trans[8]) - (trans[2] * trans[7]));
			coFact[4] = ((trans[0] * trans[8]) - (trans[2] * trans[6]));
			coFact[5] = -((trans[0] * trans[7]) - (trans[1] * trans[6]));
			coFact[6] = ((trans[1] * trans[5]) - (trans[2] * trans[4]));
			coFact[7] = -((trans[0] * trans[5]) - (trans[2] * trans[3]));
			coFact[8] = ((trans[0] * trans[4]) - (trans[1] * trans[3]));

			//update the coFact matrix to the inverse matrix by dividing each element by the determinant of the original
			for (int i = 0; i < 9; i++)
				coFact[i] = coFact[i] / detOrig;

			//Update atomic coordinates
			for (int i = 0; i < atomCoords.size(); i++)
			{
				long double holdA = atomCoords[i].a, holdB = atomCoords[i].b, holdC = atomCoords[i].c;
				atomCoords[i].a = (coFact[0] * holdA) + (coFact[1] * holdB) + (coFact[2] * holdC);
				atomCoords[i].b = (coFact[3] * holdA) + (coFact[4] * holdB) + (coFact[5] * holdC);
				atomCoords[i].c = (coFact[6] * holdA) + (coFact[7] * holdB) + (coFact[8] * holdC);
			}

			//Update tags
			directTag = true;
			cartesianTag = false;
			fetchSupercellVolume();
		}
}

//Converts one atom to direct
void Poscar::convertToDirect(Coords &crds)
{
	///apply universal scale factor to original supercell vectors, find the determinant, and transpose of the matrix made by them
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = superCellVectorA[i] * universalScaleFactor;
		superCellVectorB[i] = superCellVectorB[i] * universalScaleFactor;
		superCellVectorC[i] = superCellVectorC[i] * universalScaleFactor;
	}
	for (int i = 0; i < atomCoords.size(); i++)
	{
		atomCoords[i].a = atomCoords[i].a * universalScaleFactor;
		atomCoords[i].b = atomCoords[i].b * universalScaleFactor;
		atomCoords[i].c = atomCoords[i].c * universalScaleFactor;
	}
	universalScaleFactor = 1.0;

	//Determinant:
	long double detOrig = (superCellVectorA[0] * ((superCellVectorB[1] * superCellVectorC[2]) - (superCellVectorB[2] * superCellVectorC[1])))
		- (superCellVectorB[0] * ((superCellVectorA[1] * superCellVectorC[2]) - (superCellVectorA[2] * superCellVectorC[1])))
		+ (superCellVectorC[0] * ((superCellVectorB[2] * superCellVectorA[1]) - (superCellVectorB[1] * superCellVectorA[2])));
	if (detOrig == 0)
	{
		std::cout << "converToCartesian(): error - determinant was found to be zero\n";
		return;
	}

	//Transpose of the matrix used to convert from direct to cartesian turns out to just be the matrix of vectors as written in the head of the input file:
	long double trans[9];
	for (int i = 0; i < 3; i++)
		trans[i] = superCellVectorA[i];
	for (int i = 3; i < 6; i++)
		trans[i] = superCellVectorB[i - 3];
	for (int i = 6; i < 9; i++)
		trans[i] = superCellVectorC[i - 6];

	//Fill the cofactor matrix with the determinants of the transposed matrix, flipping the signs when necessary
	long double coFact[9];
	coFact[0] = ((trans[4] * trans[8]) - (trans[5] - trans[7]));
	coFact[1] = -((trans[3] * trans[8]) - (trans[5] - trans[6]));
	coFact[2] = ((trans[3] * trans[7]) - (trans[4] - trans[6]));
	coFact[3] = -((trans[1] * trans[8]) - (trans[2] * trans[7]));
	coFact[4] = ((trans[0] * trans[8]) - (trans[2] * trans[6]));
	coFact[5] = -((trans[0] * trans[7]) - (trans[1] * trans[6]));
	coFact[6] = ((trans[1] * trans[5]) - (trans[2] * trans[4]));
	coFact[7] = -((trans[0] * trans[5]) - (trans[2] * trans[3]));
	coFact[8] = ((trans[0] * trans[4]) - (trans[1] * trans[3]));

	//update the coFact matrix to the inverse matrix by dividing each element by the determinant of the original
	for (int i = 0; i < 9; i++)
		coFact[i] = coFact[i] / detOrig;

	//Update atomic coordinates
	long double holdA = crds.a, holdB = crds.b, holdC = crds.c;
	crds.a = (coFact[0] * holdA) + (coFact[1] * holdB) + (coFact[2] * holdC);
	crds.b = (coFact[3] * holdA) + (coFact[4] * holdB) + (coFact[5] * holdC);
	crds.c = (coFact[6] * holdA) + (coFact[7] * holdB) + (coFact[8] * holdC);


}

void Poscar::convertToDirect(std::vector <Coords> &crds)
{
	///apply universal scale factor to original supercell vectors, find the determinant, and transpose of the matrix made by them
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = superCellVectorA[i] * universalScaleFactor;
		superCellVectorB[i] = superCellVectorB[i] * universalScaleFactor;
		superCellVectorC[i] = superCellVectorC[i] * universalScaleFactor;
	}
	for (int i = 0; i < atomCoords.size(); i++)
	{
		atomCoords[i].a = atomCoords[i].a * universalScaleFactor;
		atomCoords[i].b = atomCoords[i].b * universalScaleFactor;
		atomCoords[i].c = atomCoords[i].c * universalScaleFactor;
	}
	universalScaleFactor = 1.0;

	//Determinant:
	long double detOrig = (superCellVectorA[0] * ((superCellVectorB[1] * superCellVectorC[2]) - (superCellVectorB[2] * superCellVectorC[1])))
		- (superCellVectorB[0] * ((superCellVectorA[1] * superCellVectorC[2]) - (superCellVectorA[2] * superCellVectorC[1])))
		+ (superCellVectorC[0] * ((superCellVectorB[2] * superCellVectorA[1]) - (superCellVectorB[1] * superCellVectorA[2])));
	if (detOrig == 0)
	{
		std::cout << "converToCartesian(): error - determinant was found to be zero\n";
		return;
	}

	//Transpose of the matrix used to convert from direct to cartesian turns out to just be the matrix of vectors as written in the head of the input file:
	long double trans[9];
	for (int i = 0; i < 3; i++)
		trans[i] = superCellVectorA[i];
	for (int i = 3; i < 6; i++)
		trans[i] = superCellVectorB[i - 3];
	for (int i = 6; i < 9; i++)
		trans[i] = superCellVectorC[i - 6];

	//Fill the cofactor matrix with the determinants of the transposed matrix, flipping the signs when necessary
	long double coFact[9];
	coFact[0] = ((trans[4] * trans[8]) - (trans[5] - trans[7]));
	coFact[1] = -((trans[3] * trans[8]) - (trans[5] - trans[6]));
	coFact[2] = ((trans[3] * trans[7]) - (trans[4] - trans[6]));
	coFact[3] = -((trans[1] * trans[8]) - (trans[2] * trans[7]));
	coFact[4] = ((trans[0] * trans[8]) - (trans[2] * trans[6]));
	coFact[5] = -((trans[0] * trans[7]) - (trans[1] * trans[6]));
	coFact[6] = ((trans[1] * trans[5]) - (trans[2] * trans[4]));
	coFact[7] = -((trans[0] * trans[5]) - (trans[2] * trans[3]));
	coFact[8] = ((trans[0] * trans[4]) - (trans[1] * trans[3]));

	//update the coFact matrix to the inverse matrix by dividing each element by the determinant of the original
	for (int i = 0; i < 9; i++)
		coFact[i] = coFact[i] / detOrig;

	//Update atomic coordinates
	for (int i = 0; i < crds.size(); i++)
	{
		long double holdA = crds[i].a, holdB = crds[i].b, holdC = crds[i].c;
		crds[i].a = (coFact[0] * holdA) + (coFact[1] * holdB) + (coFact[2] * holdC);
		crds[i].b = (coFact[3] * holdA) + (coFact[4] * holdB) + (coFact[5] * holdC);
		crds[i].c = (coFact[6] * holdA) + (coFact[7] * holdB) + (coFact[8] * holdC);
	}

}


//Converts all atoms to cartesian.  Could do just one atom, but I cant see a reason that you'd ever want to mix the two coordinate systems in a single file
void Poscar::convertToCartesian()
{
	if (cartesianTag) ///leaving this here just in case you'd want to do something special if file is already in cartesian
	{
		return;
	}
	else 
		if (directTag)
		{
		///conversion math
			for (int i = 0; i < atomCoords.size(); i++)
			{
				long double holdA = atomCoords[i].a, holdB = atomCoords[i].b, holdC = atomCoords[i].c;
				atomCoords[i].a = universalScaleFactor * ((superCellVectorA[0] * holdA) + (superCellVectorB[0] * holdB) + (superCellVectorC[0] * holdC));
				atomCoords[i].b = universalScaleFactor * ((superCellVectorA[1] * holdA) + (superCellVectorB[1] * holdB) + (superCellVectorC[1] * holdC));
				atomCoords[i].c = universalScaleFactor * ((superCellVectorA[2] * holdA) + (superCellVectorB[2] * holdB) + (superCellVectorC[2] * holdC));
			}

			//Updating tags
			for (int i = 0; i < 3; i++)
			{
				superCellVectorA[i] = universalScaleFactor * superCellVectorA[i];
				superCellVectorB[i] = universalScaleFactor * superCellVectorB[i];
				superCellVectorC[i] = universalScaleFactor * superCellVectorC[i];
			}

			universalScaleFactor = 1.0;
			directTag = false;
			cartesianTag = true;
			fetchSupercellVolume();
		}
}

//Converts one atom to cartesian
void Poscar::convertToCartesian(Coords &crds)
{
	///conversion math
	long double holdA = crds.a, holdB = crds.b, holdC = crds.c;
	crds.a = universalScaleFactor * ((superCellVectorA[0] * holdA) + (superCellVectorB[0] * holdB) + (superCellVectorC[0] * holdC));
	crds.b = universalScaleFactor * ((superCellVectorA[1] * holdA) + (superCellVectorB[1] * holdB) + (superCellVectorC[1] * holdC));
	crds.c = universalScaleFactor * ((superCellVectorA[2] * holdA) + (superCellVectorB[2] * holdB) + (superCellVectorC[2] * holdC));


	//Updating tags
	for (int i = 0; i < 3; i++)
	{
		superCellVectorA[i] = universalScaleFactor * superCellVectorA[i];
		superCellVectorB[i] = universalScaleFactor * superCellVectorB[i];
		superCellVectorC[i] = universalScaleFactor * superCellVectorC[i];
	}

	universalScaleFactor = 1.0;
}

//Converts all atoms to cartesian.  
void Poscar::convertToCartesian(std::vector <Coords> &crds)
{
			///conversion math
			for (int i = 0; i < crds.size(); i++)
			{
				long double holdA = crds[i].a, holdB = crds[i].b, holdC = crds[i].c;
				crds[i].a = universalScaleFactor * ((superCellVectorA[0] * holdA) + (superCellVectorB[0] * holdB) + (superCellVectorC[0] * holdC));
				crds[i].b = universalScaleFactor * ((superCellVectorA[1] * holdA) + (superCellVectorB[1] * holdB) + (superCellVectorC[1] * holdC));
				crds[i].c = universalScaleFactor * ((superCellVectorA[2] * holdA) + (superCellVectorB[2] * holdB) + (superCellVectorC[2] * holdC));
			}

			//Updating tags
			for (int i = 0; i < 3; i++)
			{
				superCellVectorA[i] = universalScaleFactor * superCellVectorA[i];
				superCellVectorB[i] = universalScaleFactor * superCellVectorB[i];
				superCellVectorC[i] = universalScaleFactor * superCellVectorC[i];
			}

			universalScaleFactor = 1.0;
}

//Removes atoms that have their 'extraInfo' marked as 'tag'.
void Poscar::removeTaggedAtoms(std::string tag)
{
	std::vector <Coords> tmp;

	for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo != tag)
					tmp.push_back(atomCoords[i]);

	atomCoords = tmp;
	updateAtomCounts(tmp);
}

//Removes atoms that have the same R3 coordinates.  does NOT take any other parts of the Coords class into account (id, extraInfo, etc.)
void Poscar::removeDuplicates()
{
	std::vector <Coords> tmp;

	for (int i = 0; i < atomCoords.size(); i++)
		for (int j = 0; j < atomCoords.size(); j++)
			if (i != j)
				if (atomCoords[i] == atomCoords[j] && atomCoords[i].extraInfo != "cp")
					atomCoords[j].extraInfo = "cp";
				
	for (int i = 0; i < atomCoords.size(); i++)
		if (atomCoords[i].extraInfo != "cp")
			tmp.push_back(atomCoords[i]);

	atomCoords = tmp;
	updateAtomCounts(tmp);
}

void Poscar::removeDuplicates(double tol)
{
	std::vector <Coords> tmp;

	for (int i = 0; i < atomCoords.size(); i++)
		for (int j = 0; j < atomCoords.size(); j++)
			if (i != j)
				if ((dist_(atomCoords[i], atomCoords[j]) < tol) && (atomCoords[i].extraInfo != "cp"))
					atomCoords[j].extraInfo = "cp";

	for (int i = 0; i < atomCoords.size(); i++)
		if (atomCoords[i].extraInfo != "cp")
			tmp.push_back(atomCoords[i]);

	atomCoords = tmp;
	updateAtomCounts(tmp);
}

void Poscar::removeDuplicates(std::vector <Coords> &atoms )
{
	std::vector <Coords> tmp;

	for (int i = 0; i < atoms.size(); i++)
		for (int j = 0; j < atoms.size(); j++)
			if (i != j)
				if (atoms[i] == atoms[j] && atoms[i].extraInfo != "cp")
					atoms[j].extraInfo = "cp";

	for (int i = 0; i < atoms.size(); i++)
		if (atoms[i].extraInfo != "cp")
			tmp.push_back(atoms[i]);

	atoms = tmp;
}

void Poscar::removeDuplicates(std::string instruct)
{
	std::vector <Coords> tmp;

	if (instruct == "coords")
	{
		for (int i = 0; i < atomCoords.size(); i++)
			for (int j = 0; j < atomCoords.size(); j++)
				if (i != j)
					if (atomCoords[i] == atomCoords[j] && atomCoords[i].extraInfo != "cp")
						atomCoords[j].extraInfo = "cp";

		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo != "cp")
				tmp.push_back(atomCoords[i]);
	}

	if (instruct == "id")
	{
		for (int i = 0; i < atomCoords.size(); i++)
			for (int j = 0; j < atomCoords.size(); j++)
				if (i != j)
					if (atomCoords[i].id == atomCoords[j].id && atomCoords[i].extraInfo != "cp")
						atomCoords[j].extraInfo = "cp";

		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo != "cp")
				tmp.push_back(atomCoords[i]);
	}

	if (instruct == "unoriginal")
	{
		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo != "original")
				atomCoords[i].extraInfo = "cp";

		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo != "cp")
				tmp.push_back(atomCoords[i]);
	}

	if (instruct == "edge")
	{
		convertToDirect();
		for (int i = 0; i < atomCoords.size(); i++)
		{
			if (atomCoords[i].a > 0.9999 || atomCoords[i].b > 0.9999 || atomCoords[i].c > 0.9999)
				atomCoords[i].extraInfo = "rm";
			if (atomCoords[i].a < -0.0001 || atomCoords[i].b < -0.0001 || atomCoords[i].c < -0.0001)
				atomCoords[i].extraInfo = "rm";
		}

		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].extraInfo != "rm")
				tmp.push_back(atomCoords[i]);
	}

	atomCoords = tmp;
	updateAtomCounts(tmp);
}

void Poscar::updateAtomTypes()
{
	bool ok = true;
	std::vector <std::string> tmp, tmp2;

	for (int i = 0; i < atomCoords.size(); i++)
		if (atomCoords[i].atomType != "undf" && atomCoords[i].atomType != "END" && (atomCoords[i].atomType[0] != 'u' && atomCoords[i].atomType[1] != 'n'))
			tmp.push_back(atomCoords[i].atomType);

	std::sort(tmp.begin(), tmp.end(), stringSort_);
	tmp.push_back("END");

	for (int i = 0; i < tmp.size(); i++)
		if (tmp[i] != "END")
			if (tmp[i] != tmp[i+1])
				tmp2.push_back(tmp[i]);

	atomTypes = tmp2;
}

void Poscar::updateAtomCounts()
{
	std::vector <int> tmp;
	for (int i = 0; i < atomTypes.size(); i++)
		tmp.push_back(0);


	for (int i = 0; i < tmp.size(); i++)
		for (int j = 0; j < atomCoords.size(); j++)
			if (atomCoords[j].atomType == atomTypes[i])
				tmp[i] ++;

	atomTypeNums = tmp;
}

void Poscar::updateAtomCounts(std::vector <Coords> crds)
{
	//Update params
	for (int i = 0; i < atomTypeNums.size(); i++)
		atomTypeNums[i] = 0;

	for (int i = 0; i < atomTypeNums.size(); i++)
		for (int j = 0; j < crds.size(); j++)
			if (crds[j].atomType == atomTypes[i])
				atomTypeNums[i] ++;
}

void Poscar::updateAtomCoords()
{
	if (directTag)
		convertToDirect();
	if (cartesianTag)
		convertToCartesian();
}

void Poscar::assertAtomOrder(std::vector <std::string> order)
{
	std::vector <std::string> newOrder;
	for (int i = 0; i < order.size(); i++){
		bool counted = false;
		for(int j = 0; j < newOrder.size(); j++)
			if(order[i] == newOrder[j])
				counted = true;
		if(!counted)
			newOrder.push_back(order[i]);
	}

	std::vector <Coords> tmp;
	for (int i = 0; i < newOrder.size(); i++)
		for (int j = 0; j < atomCoords.size(); j++)
			if(newOrder[i] == atomCoords[j].atomType && atomCoords[j].extraInfo != "counted"){
				tmp.push_back(atomCoords[j]);
				atomCoords[j].extraInfo = "counted";
			}
		
	if(tmp.size() != atomCoords.size()){
		std::cout << "WARNING: AssertAtomOrder: New atom coords size does not match old atom coords size!\n";
		std::cout << "I REFUSE TO CONTINUE WITH THIS SICK JOB...BYE!!!";
	}	
	else{
		atomTypes = newOrder;
		atomCoords = tmp;
		updateAtomCounts();
	}
	return;
}

void Poscar::updateAll()
{
	updateAtomTypes();
	updateAtomCounts();
	updateAtomCoords();
}

bool atomCoordsSortRules(Coords a_, Coords b_)
{
	return a_.id < b_.id;
}

void Poscar::copyFormatting(const Poscar toCopy)
{
	std::vector <std::string> orderForStrings;
	std::vector <std::vector <Coords> > atomsByType;

	//Get the order that the atom types should be in
	for (int i = 0; i < toCopy.atomTypes.size(); i++)
		orderForStrings.push_back(toCopy.atomTypes[i]);

	std::vector <std::string> extraTypes;
	//Re-Order the current atom types, if theyre all there.  Add whatever extras there might be at the end
	for (int i = 0; i < atomTypes.size(); i++)
	{
		bool currentIsThere = false; ///assume current considered atom does not exist in the new coords
		for (int j = 0; j < toCopy.atomTypes.size(); j++)
			if (atomTypes[i] == toCopy.atomTypes[j])
				currentIsThere = true; ///if it does exist, say so
		if (!currentIsThere)
			extraTypes.push_back(atomTypes[i]); ///but if it dosent, push the atom type into the extras
	}

	//Update current atom types, adding extras to the end
	atomTypes = toCopy.atomTypes;
	for (int i = 0; i < extraTypes.size(); i++)
		atomTypes.push_back(extraTypes[i]);

	std::vector <Coords> tmpCoords;
	std::vector <Coords> extraCoords;
	//Sort the atom coords by ID, to match the original
	///Sort the current atoms by this order, adding extras to the end
	for (int i = 0; i < toCopy.atomCoords.size(); i++)
	{
		bool idIsThere = false;
		for (int j = 0; j < atomCoords.size(); j++)
			if (toCopy.atomCoords[i].id == atomCoords[j].id)
			{
				idIsThere = true;
				tmpCoords.push_back(atomCoords[j]);
			}
		if (!idIsThere)
			extraCoords.push_back(atomCoords[i]);
	}

	//Update current atom counts, adding extras to the end
	atomCoords = tmpCoords;
	for (int i = 0; i < extraCoords.size(); i++)
		atomCoords.push_back(extraCoords[i]);

	//Finially, update all of the atom type counts and hope it dosent reverse everything that I've just done.  Also switch coord systems
	updateAtomCounts();
	if (toCopy.directTag)
		convertToDirect();
	else
		if (toCopy.cartesianTag)
			convertToCartesian();
}

void Poscar::print()
{
	//Print original input file
	using namespace std;
	cout << "Printing input file..." << endl << endl << endl;

	//Open file
	ifstream infile(infilePath.c_str());
	if (infile.fail())
		cout << "print(): File could not be found...expect errors\n";

	//Print line by line
	string s;
	while (!infile.eof())
	{
		getline(infile, s);
		cout << s << endl;
	}

	infile.close();

	//Print class contents (mostly for debugging)
	cout << endl << endl << "Printing class interpretation of input file..." << endl << endl << endl;
	cout.setf(ios::fixed);
	cout.setf(ios::showpoint);

	cout << fileTitle << endl;
	cout << universalScaleFactor << endl;
	cout << left << setfill('0') << setprecision(10) << superCellVectorA[0] << "\t" << superCellVectorA[1] << "\t" << superCellVectorA[2] << endl;
	cout << left << setfill('0') << setprecision(10) << superCellVectorB[0] << "\t" << superCellVectorB[1] << "\t" << superCellVectorB[2] << endl;
	cout << left << setfill('0') << setprecision(10) << superCellVectorC[0] << "\t" << superCellVectorC[1] << "\t" << superCellVectorC[2] << endl;
	for (int i = 0; i < atomTypes.size(); i++)
		cout << atomTypes[i] << "   ";
	cout << endl;
	for (int i = 0; i < atomTypeNums.size(); i++)
		cout << atomTypeNums[i] << "   ";
	cout << endl;
	if (selectiveDynamicsTag)
		cout << "Sel" << endl;
	if (directTag)
		cout << "Direct   ";
	if (cartesianTag)
		cout << "Cartesian   ";
	cout << endl;
	for (int j = 0; j < atomTypes.size(); j++)
		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].atomType == atomTypes[j])
				cout << left << setfill('0') << setprecision(10) << atomCoords[i].a << "\t" << atomCoords[i].b << "\t" << atomCoords[i].c << "\t" << atomCoords[i].flags << endl;
	cout << "FILE END.  THERE SHOULD BE NO SPACE BETWEEN THIS MESSAGE AND THE FINAL LINE OF THE FILE." << endl;
}

void Poscar::print(std::string headOrFull, std::string normalOrExtra)
{
	using namespace std;
	//Always print head of file
	cout << endl << endl << "Printing class interpretation of input file..." << endl << endl << endl;
	cout.setf(ios::fixed);
	cout.setf(ios::showpoint);
	
	if (normalOrExtra != "extra" && normalOrExtra != "Extra")
	{
		cout << fileTitle << endl;
		cout << universalScaleFactor << endl;
		cout << left << setfill('0') << setprecision(10) << superCellVectorA[0] << "\t" << superCellVectorA[1] << "\t" << superCellVectorA[2] << endl;
		cout << left << setfill('0') << setprecision(10) << superCellVectorB[0] << "\t" << superCellVectorB[1] << "\t" << superCellVectorB[2] << endl;
		cout << left << setfill('0') << setprecision(10) << superCellVectorC[0] << "\t" << superCellVectorC[1] << "\t" << superCellVectorC[2] << endl;
		for (int i = 0; i < atomTypes.size(); i++)
			cout << atomTypes[i] << "   ";
		cout << endl;
		for (int i = 0; i < atomTypeNums.size(); i++)
			cout << atomTypeNums[i] << "   ";
		cout << endl;
		if (selectiveDynamicsTag)
			cout << "Sel" << endl;
		if (directTag)
			cout << "Direct   ";
		if (cartesianTag)
			cout << "Cartesian   ";
		cout << endl;
	}

	int linenum = 1;
	if (normalOrExtra == "extra" || normalOrExtra == "Extra")
	{
		cout << "Linenum " << linenum << "\t" << fileTitle << endl;
		linenum ++;
		cout << "Linenum " << linenum << "\t" << universalScaleFactor << endl;
		linenum++;
		cout << "Linenum " << linenum << "\t" << left << setfill('0') << setprecision(10) << superCellVectorA[0] << "\t" << superCellVectorA[1] << "\t" << superCellVectorA[2] << endl;
		linenum++;
		cout << "Linenum " << linenum << "\t" << left << setfill('0') << setprecision(10) << superCellVectorB[0] << "\t" << superCellVectorB[1] << "\t" << superCellVectorB[2] << endl;
		linenum++;
		cout << "Linenum " << linenum << "\t" << left << setfill('0') << setprecision(10) << superCellVectorC[0] << "\t" << superCellVectorC[1] << "\t" << superCellVectorC[2] << endl;
		linenum++;
		cout << "Linenum " << linenum << "\t";
		linenum++;
		for (int i = 0; i < atomTypes.size(); i++)
			cout << atomTypes[i] << "   ";
		cout << endl;
		cout << "Linenum " << linenum << "\t";
		for (int i = 0; i < atomTypeNums.size(); i++)
			cout << atomTypeNums[i] << "   ";
		cout << endl;
		if (selectiveDynamicsTag)
		{
			linenum++;
			cout << "Linenum " << linenum << "\tSel" << endl;
		}
		linenum++;
		cout << "Linenum " << linenum << "\t";
		if (directTag)
			cout << "Direct   ";
		if (cartesianTag)
			cout << "Cartesian   ";
		cout << endl;
	}

	//if user wants the rest of file, print that.  check for extra options as well
	if (headOrFull == "Full" || headOrFull == "full")
	{
		for (int j = 0; j < atomTypes.size(); j++)
		{
			for (int i = 0; i < atomCoords.size(); i++)
				if (atomCoords[i].atomType == atomTypes[j])
				{
					if (normalOrExtra == "extra" || normalOrExtra == "Extra")
					{
						linenum++;
						cout << "Linenum " << linenum << "\tAtomnum" << i + 1 << "\t";
					}
					cout << left << setfill('0') << setprecision(10) << atomCoords[i].a << "\t" << atomCoords[i].b << "\t" << atomCoords[i].c << "\t" << atomCoords[i].flags;
					if (normalOrExtra == "extra" || normalOrExtra == "Extra")
						cout << "\tAtom Type: " << atomCoords[i].atomType << "   extra:  " << atomCoords[i].extraInfo << endl;
					cout << endl;
				}
		}
		cout << "FILE END.  THERE SHOULD BE NO SPACE BETWEEN THIS MESSAGE AND THE FINAL LINE OF THE FILE." << endl;
	}

}

void Poscar::write()
{
	//Open file
	using namespace std;
	ofstream outfile(defaultWritePath.c_str());


	//Write file.  literally the exact same as print, except 'cout' has been replaced with 'outfile'
	outfile << fileTitle << endl;

	outfile << fixed;
	outfile << setprecision(3);

	outfile << left << universalScaleFactor << endl;
	outfile << setprecision(10);
	outfile << left << superCellVectorA[0] << "   " << superCellVectorA[1] << "   " << superCellVectorA[2] << endl;
	outfile << left << superCellVectorB[0] << "   " << superCellVectorB[1] << "   " << superCellVectorB[2] << endl;
	outfile << left << superCellVectorC[0] << "   " << superCellVectorC[1] << "   " << superCellVectorC[2] << endl;
	for (int i = 0; i < atomTypes.size(); i++)
		outfile << atomTypes[i] << "   ";
	outfile << endl;
	for (int i = 0; i < atomTypeNums.size(); i++)
		outfile << atomTypeNums[i] << "   ";
	outfile << endl;
	if (selectiveDynamicsTag)
		outfile << "Sel" << endl;
	if (directTag)
		outfile << "Direct   ";
	if (cartesianTag)
		outfile << "Cartesian   ";
	outfile << endl;
	//Print atoms in order, just in case they got mixed up somehow
	for (int j = 0; j < atomTypes.size(); j++)
		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].atomType == atomTypes[j])
				outfile << left << setfill('0') << setprecision(10) << atomCoords[i].a << "   " << atomCoords[i].b << "   "
				<< atomCoords[i].c << "   " << atomCoords[i].flags << endl;
	
	outfile.close();
}

void Poscar::write(std::string userPath)
{
	//Open file
	using namespace std;
	ofstream outfile(userPath.c_str());

	//Write file.  literally the exact same as print, except 'cout' has been replaced with 'outfile'
	outfile << fileTitle << endl;

	outfile << fixed;
	outfile << setprecision(3);

	outfile << left << universalScaleFactor << endl;
	outfile << setprecision(10);
	outfile << left << superCellVectorA[0] << "   " << superCellVectorA[1] << "   " << superCellVectorA[2] << endl;
	outfile << left << superCellVectorB[0] << "   " << superCellVectorB[1] << "   " << superCellVectorB[2] << endl;
	outfile << left << superCellVectorC[0] << "   " << superCellVectorC[1] << "   " << superCellVectorC[2] << endl;
	for (int i = 0; i < atomTypes.size(); i++)
		outfile << atomTypes[i] << "   ";
	outfile << endl;
	for (int i = 0; i < atomTypeNums.size(); i++)
		outfile << atomTypeNums[i] << "   ";
	outfile << endl;
	if (selectiveDynamicsTag)
		outfile << "Sel" << endl;
	if (directTag)
		outfile << "Direct   ";
	if (cartesianTag)
		outfile << "Cartesian   ";
	outfile << endl;
	//Print atoms in order, just in case they got mixed up somehow
	for (int j = 0; j < atomTypes.size(); j++)
		for (int i = 0; i < atomCoords.size(); i++)
			if (atomCoords[i].atomType == atomTypes[j])
				outfile << left << setfill('0') << setprecision(10) << atomCoords[i].a << "\t" << atomCoords[i].b << "\t" << atomCoords[i].c << "\t" << atomCoords[i].flags << endl;

	outfile.close();
}

void Poscar::clearEmptyValues()
{
	for (int i = 0; i < atomTypes.size(); i++)
		if (atomTypes[i] == "" || atomTypes[i] == " ")
			atomTypes[i] = "del";
	std::vector <std::string> tmp;

	for (int i = 0; i < atomTypes.size(); i++)
		if (atomTypes[i] != "del")
			tmp.push_back(atomTypes[i].c_str());

	atomTypes = tmp;
}
