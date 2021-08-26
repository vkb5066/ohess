#Deals with computing Cij values from OHESS strains
#OHESS formalism given in 'High-efficiency calculation of elastic constants enhanced by the optimized
#                          strain-matrix sets, Zhong-Li Liu (2020)'

#Works in reference to the definitions and syntax from my C++ deform code
#Syntax from C++ code:

#//Database Keys:
#/*
#* 0 = Cubic
#* 1 = Hexagonal
#* 2 = Rhomb I
#* 3 = Rhomb II
#* 4 = Tetrag I
#* 5 = Tetrag II
#* 6 = Ortho
#* 7 = Mono
#* 8 = Tri
#*/

#extern const double OHESS_DB[9][6][6] =
#{	//---CUBIC--- (1 nonzero)
#	{
#		{1.0, 0.0, 0.0, 1.0, 0.0, 0.0}, //c11, c12, c44
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
#	},
#	//---HEXAG--- (2 nonzeros)
#	{
#		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c11, c12, c13
#		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c44
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
#	},
#	//---RHOMB I--- (2 nonzeros)
#	{
#		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c1i, i = {1-4}
#		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c44
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
#	},
#	//---RHOMB II--- (2 nonzeros)
#	{
#		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c1i, i = {1-5}
#		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c44
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
#	},
#	//---TETRAG I--- (2 nonzeros)
#	{
#		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c11, c12, c13
#		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c44
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
#	},
#	//---TETRAG II--- (2 nonzeros)
#	{
#		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c1i, i = {1-3, 6}
#		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c44
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
#	},
#	//---ORTHO--- (3 nonzeros)
#	{
#		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c11, c12, c13
#		{0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, //c22, c23
#		{0.0, 0.0, 1.0, 1.0, 1.0, 1.0}, //cii, i = {3-6}
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
#	},
#	//---MONO--- (4 nonzeros)
#	{
#		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c1i, i = {1-3, 6}
#		{0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, //c22, c23, c26
#		{0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, //c33, c36, c44, c45
#		{0.0, 0.0, 0.0, 0.0, 1.0, 1.0}, //c55, c66
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
#	},
#	//---TRICLIN--- (6 nonzeros)
#	{
#		{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}, //c1i, i = {1-6}
#		{0.0, 1.0, 0.0, 0.0, 0.0, 0.0}, //c2i, i = {2-6}
#		{0.0, 0.0, 1.0, 0.0, 0.0, 0.0}, //c3i, i = {3-6}
#		{0.0, 0.0, 0.0, 1.0, 0.0, 0.0}, //c44, c45, c46
#		{0.0, 0.0, 0.0, 0.0, 1.0, 0.0}, //c55, c56
#		{0.0, 0.0, 0.0, 0.0, 0.0, 1.0}  //c66
#	}
#};

from scipy.stats import linregress
from copy import deepcopy

ZERO_VOIGT_MAP = {0: "xx", 1: "yy", 2: "zz", 3: "yz", 4: "zx", 5: "xy"}
NORM_VOIGT_MAP = {1: "xx", 2: "yy", 3: "zz", 4: "yz", 5: "zx", 6: "xy"}

def KbToGPa(inKb, reverse=False):
    return -0.1 * inKb if not reverse else -10.0 * inKb

def AlmEqu(a, b, tol=1e-7):
    return (abs(a - b) < tol)

#Converts a delta (defined in ref paper) to the fractional strain
#OHESS seems to be much simpler in its definitions of delta than the standard E vs strain methods
#I'm keeping this here in case I am wrong
def DeltaToEpsilon(delta):
    return delta

#Returns the {(i, j)} that can be calculated from a given strain set in a given lattice system
#!!! Everything is zero-indexed !!! i.e. voigt's C23 is this function's C12
#Only returns one of the ij, symmetry is implicit
def GetAllowedCijs(latSysId, strainMatrixId):
    ##Cubic
    if(latSysId == 0):
        if(strainMatrixId == 0):
            return [[0, 0], [0, 1], [3, 3]]
    ##Hexag
    if(latSysId == 1):
        if(strainMatrixId == 0):
            return [[0, i] for i in range(0, 3)]
        if(strainMatrixId == 1):
            return [[2, 2], [3, 3]]
    ##Rhomb I
    if(latSysId == 2):
        if(strainMatrixId == 0):
            return [[0, i] for i in range(0, 4)]
        if(strainMatrixId == 1):
            return [[2, 2], [3, 3]]
    ##Rhomb II
    if(latSysId == 3):
        if(strainMatrixId == 0):
            return [[0, i] for i in range(0, 5)]
        if(strainMatrixId == 1):
            return [[2, 2], [3, 3]]
    ##Tetrag I
    if(latSysId == 4):
        if(strainMatrixId == 0):
            return [[0, i] for i in range(0, 3)]
        if(strainMatrixId == 1):
            return [[2, 2], [3, 3]]
    ##Tetrag II
    if(latSysId == 5):
        if(strainMatrixId == 0):
            return [[0, i] for i in range(0, 3)] + [[0, 5]]
        if(strainMatrixId == 1):
            return [[2, 2], [3, 3]]
    ##Ortho
    if(latSysId == 6):
        if(strainMatrixId == 0):
            return [[0, i] for i in range(0, 3)]
        if(strainMatrixId == 1):
            return [[1, 1], [1, 2]]
        if(strainMatrixId == 2):
            return [[i, i] for i in range(2, 6)]
    ##Mono
    if(latSysId == 7):
        if(strainMatrixId == 0):
            return [[0, i] for i in range(0, 3)] + [[0, 5]]
        if(strainMatrixId == 1):
            return [[1, 1], [1, 2], [1, 5]]
        if(strainMatrixId == 2):
            return [[2, 2], [2, 5], [3, 3], [3, 4]]
        if(strainMatrixId == 3):
            return [[4, 4], [5, 5]]
    ##Tricl
    if(latSysId == 8):
        return [[strainMatrixId, i] for i in range(strainMatrixId, 6)]

    ##Otherwise, return empty set
    return []

#From 'Physical Properties of Crystals', Nye (1985), pgs 140 - 141 (mostly) and
#'Irreducible matrix resolution of the elasticity tensor for symmetry systems', Yakov Itin, pg 19
#Returns the UPPER TRIANGLE portion's symmetric Cijs considering crystal class !!! (zero-indexed) !!!
#General, but not complete rules.  Only works if following the OHESS method's allowed Cijs
#Otherwise, using this will set components of the matrix to zero.  So be careful of that
def GetEquivCijs(latSysId, c):
    ##Cubic
    if(latSysId == 0):
        c[1][1], c[2][2] = c[0][0], c[0][0] ###C11 = C22 = C33
        c[0][2], c[1][2] = c[0][1], c[0][1] ###C12 = C13 = C23
        c[4][4], c[5][5] = c[3][3], c[3][3] ###C44 = C55 = C66
    ##Hexagonal
    if(latSysId == 1):
        c[1][1] = c[0][0] ###C11 = C22
        c[1][2] = c[0][2] ###C13 = C23
        c[4][4] = c[3][3] ###C44 = C55
        c[5][5] = 0.5*(c[0][0] - c[0][1]) ###C66 = 1/2(C11 - C12)
    ##Rhomb I, II (aka Trigonal I, II)
    if(latSysId == 2 or latSysId == 3):
        c[1][1] = c[0][0] ###C11 = C22
        c[1][2] = c[0][2] ###C13 = C23
        c[1][3], c[4][5] = -c[0][3], c[0][3] ###C14 = -C24, C14 = C56
    ##Rhomb II (extra)
    if(latSysId == 3):
        c[1][4], c[3][5] = -c[0][4], -c[0][4] ###C15 = -C25, C15 = -C46
    ##Tetrag I, II
    if(latSysId == 4 or latSysId == 5):
        c[1][1] = c[0][0] ###C11 = C22
        c[4][4] = c[3][3] ###C44 = C55
        c[1][2] = c[0][2] ###C13 = C23
    ##Tetrag II (extra)
    if(latSysId == 5):
        c[1][5] = -c[0][5] ###C16 = -C26
    ##Ortho
    if(latSysId == 6):
        pass ###No symmetry operations
    ##Mono
    if(latSysId == 7):
        pass ###No symmetry operations
    ##Triclin
    if(latSysId == 8):
        pass ###No symmetry operations

    return

#Returns the lattice system ID, deform matrix ID, and delta from the work directory
#Work directory is ALWAYS written as N_M_sD OR N_-1_0.00000 (for no applied strain)
#Where 0 <= N < 9 is the lattice system index,
#      0 <= M < 6 is the deform matrix index,
#      s is either the character 'n' or the character 'p' for negative or positive,
#and   D is the fractional delta written to 5 decimal points
#Returns a list as [int(lat sys id), int(def mat id), float(delta)]
def ParseWorkDir(wd):
    spl = wd.split('_')

    ##No applied strain case
    if(spl[1] == "-1"):
        return [int(spl[0]), -1, 0.00000]

    signDelta = spl[2][0] ##the character n or p
    absDelta = spl[2][1:] ##the rest of delta

    return int(spl[0]), int(spl[1]), -1*float(absDelta) if signDelta == 'n' else +1*float(absDelta)

class OhessSet:
    dir = None

    latSystem = None ##0 - 8
    data = None ##in the form {matrixId: dict("epsilons": [...], "stresses": dict("xx": [...], etc)}

    cijMatrix = None
    maxErr = None ##in GPa

    #Need to initialize before actually adding any data since data for OhessSet is written as multiple lines
    def __init__(self):
        self.dir = None
        self.latSystem = None
        self.data = dict()
        self.cijMatrix = []
        for i in range(0, 6):
            self.cijMatrix.append([])
            for j in range(0, 6):
                self.cijMatrix[i].append(0.0)
        self.maxErr = 0.0

    ##'line' is written as:
    ##dir,workDir,nrg,nConv,irredKpts,aV,bV,cV,alpha,beta,gamma,vol,sxx,syy,szz,sxy,syz,szx
    def AddData(self, line):
        spl = line.replace("\n", '').split(',')
        self.dir = spl[0]

        latSysId, defMatId, delta = ParseWorkDir(spl[1])
        self.latSystem = latSysId
        ##no, the ordering isn't a mistake.  VASP writes out the stresses in a weird way
        s0, s1, s2 = float(spl[12]), float(spl[13]), float(spl[14])
        s3, s4, s5 = float(spl[16]), float(spl[17]), float(spl[15])
        if(defMatId not in self.data.keys()):
            self.data[defMatId] = {"epsilons": [DeltaToEpsilon(delta)], "stresses": {"xx": [s0], "yy": [s1],
                                                                                     "zz": [s2], "yz": [s3],
                                                                                     "zx": [s4], "xy": [s5]}}
        else:
            self.data[defMatId]["epsilons"].append(DeltaToEpsilon(delta))
            self.data[defMatId]["stresses"]["xx"].append(s0)
            self.data[defMatId]["stresses"]["yy"].append(s1)
            self.data[defMatId]["stresses"]["zz"].append(s2)
            self.data[defMatId]["stresses"]["yz"].append(s3)
            self.data[defMatId]["stresses"]["zx"].append(s4)
            self.data[defMatId]["stresses"]["xy"].append(s5)

        return

    #Applies the zero strain data to all deform matrix sets
    #Removes the -1 matrix id from data (since it corresponds to the initial system)
    #Use after all data is read in or it will definetly cause errors
    def AddZeroStrainData(self):
        zeroStrain = self.data.pop(-1)
        for ke in self.data.keys():
            self.data[ke]["epsilons"].append(zeroStrain["epsilons"][0]) ##[0] since there is only 1 entry
            for voigt in ["xx", "yy", "zz", "yz", "zx", "xy"]:
                self.data[ke]["stresses"][voigt].append(zeroStrain["stresses"][voigt][0])

        return

    def ConvertToGPa(self):
        for ke in self.data.keys():
            for keStress in self.data[ke]["stresses"]:
                self.data[ke]["stresses"][keStress] = [KbToGPa(kb) for kb in \
                                                       self.data[ke]["stresses"][keStress]]

    #i, j are !!! zero indexed !!!
    def GetCij(self, i, j):
        for defMatIndex in range(0, 6):
            if([i, j] in GetAllowedCijs(latSysId=self.latSystem, strainMatrixId=defMatIndex)):
                res = linregress(x=self.data[defMatIndex]["epsilons"],
                                 y=self.data[defMatIndex]["stresses"][ZERO_VOIGT_MAP[j]])
                ##Reasonable est of error is the diff of the intercept to zero (in GPa)
                if(abs(res[1]) > self.maxErr):
                    self.maxErr = abs(res[1])

                ##Apply Cij to this objects matrix and also return the slope
                self.cijMatrix[i][j] = res[0]

                return res[0]
        return 0.0 ##If we can't compute it, it must have been zero by symmetry

    #Get the symmetry inequivelent parts of Cij matrix
    def GetSymCijMatrix(self):
        for i in range(0, 6):
            for j in range(0, 6):
                _ = self.GetCij(i, j)

        return

    #Apply symmetry conditions to get full Cij matrix (use AFTER GetSymCijMatrix())
    def GetFullCijMatrix(self):
        ##Get upper triangular portion
        GetEquivCijs(self.latSystem, self.cijMatrix)
        ##Set Cij = Cji (symmetric stress/strain tensors)
        for i in range(0, 6):
            for j in range(i, 6):
                self.cijMatrix[j][i] = self.cijMatrix[i][j]

        return

#Returns a dict of ohess objects.  Key is set to `pwd`.split(sep)[indexToMakeKey]
def GetOhessDict(infileLoc, indexToMakeKey, sep="auto"):
    ret = {}

    with open(infileLoc, 'r') as infile:
        for num, line in enumerate(infile):
            ##Skip header
            if(num == 0):
                continue
            lin = line.split(',')

            ##If sep = auto, try to figure out what the seperator should be
            if(sep == "auto" and num == 1):
                if("/" in lin[0]):
                    sep = "/"
                elif("\\" in lin[0]):
                    sep = "\\"
                else:
                    print("Could not determine seperating character :(\n")
                    exit(1)

            ##If key not in dict, initialize a new ohess set.  Otherwise add data
            key = lin[0].split(sep)[indexToMakeKey]
            if(key not in ret.keys()):
                ret[key] = OhessSet()
            ret[key].AddData(line)
        infile.close()

    #Append initial strain data now that infile is finished being read, finish initializing
    for ke in ret.keys():
        ret[ke].AddZeroStrainData()
        ret[ke].ConvertToGPa()
        ret[ke].GetSymCijMatrix()
        ret[ke].GetFullCijMatrix()
    return ret
